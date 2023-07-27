/*
 * Copyright (c) 2016-2020, NVIDIA CORPORATION. All rights reserved.
 *
 * See COPYRIGHT for license information
 */

#include "nvshmem.h"

#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <string.h>
#include <map>
#include <algorithm>
#include "nvshmem_internal.h"
#include "nvshmemx_error.h"
#include "nvshmem_nvtx.hpp"
#include "dlmalloc.h"
#include "util.h"
#include "sockets.h"
#include "nvshmemi_team.h"

/* This is a requirement imposed by DMA-BUF which only supports 32-bit registrations */
#define NVSHMEMI_DMA_BUF_MAX_LENGTH 0x100000000ULL
#define NVSHMEMI_MAX_HANDLE_LENGTH 2147483648ULL

#define IPC_CHECK(ipcFuncResult)                                   \
    if (ipcFuncResult == -1)                                       \
    {                                                              \
        fprintf(stderr, "Failure at %u %s\n", __LINE__, __FILE__); \
        exit(EXIT_FAILURE);                                        \
    }

size_t cumem_granularity;
size_t log2_cumem_granularity;
static bool nvshmemi_no_vmm_heap_initialized;
static std::map<pid_t, int>
    p2p_processes; /* Map from p2p processes to PE id - required when using VMM */
static int is_mem_handle_null(nvshmem_mem_handle_t *handle)
{
    assert(sizeof(nvshmem_mem_handle_t) % sizeof(uint64_t) == 0);

    for (size_t i = 0; i < (sizeof(nvshmem_mem_handle_t) / sizeof(uint64_t)); i++)
    {
        if (*((uint64_t *)handle + i) != (uint64_t)0)
            return 0;
    }

    return 1;
}

static int cleanup_local_handles(nvshmem_mem_handle_t *handles, nvshmemi_state_t *state)
{
    int status = 0;

    for (int i = 0; i < state->num_initialized_transports; i++)
    {
        if (state->transport_bitmap & (1 << i))
        {
            if (!is_mem_handle_null(&handles[i]))
            {
                status = state->transports[i]->host_ops.release_mem_handle(&handles[i],
                                                                           state->transports[i]);
                NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                                      "transport release memhandle failed \n");
            }
        }
    }

out:
    return status;
}

template <typename T>
int check_for_symmetry(T value)
{
    int status = 0;
    nvshmemi_state_t *state = nvshmemi_state;

    /*TODO: need to handle multi-threaded scenarios*/

    if (!nvshmemi_options.ENABLE_ERROR_CHECKS)
        return 0;

    if (state->scratch_size < sizeof(T) * state->npes)
    {
        if (state->scratch_size)
            free(state->scratch_space);

        state->scratch_space = (char *)malloc(sizeof(T) * state->npes);
        NVSHMEMI_NULL_ERROR_JMP(state->scratch_space, status, NVSHMEMX_ERROR_OUT_OF_MEMORY, out,
                                "error allocating scratch space \n");
        state->scratch_size = sizeof(T) * state->npes;
    }

    status = nvshmemi_boot_handle.allgather((void *)&value, (void *)state->scratch_space, sizeof(T),
                                            &nvshmemi_boot_handle);
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                          "allgather in symmetry check failed \n");

    for (int i = 0; i < state->npes; i++)
    {
        status = (*((T *)state->scratch_space + i) == value) ? 0 : 1;
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_SYMMETRY, out, "symmetry check failed \n");
    }

out:
    return status;
}

int mspace_track_large_chunks(mspace msp, int enable);
size_t destroy_mspace(mspace msp);

int nvshmemi_setup_memory_space(nvshmemi_state_t *state, void *heap_base, size_t size)
{
    int status = 0;
    mspace heap_mspace = 0;

    heap_mspace = create_mspace_with_base(heap_base, size, 0);
    NVSHMEMI_NULL_ERROR_JMP(heap_mspace, status, NVSHMEMX_ERROR_OUT_OF_MEMORY, out,
                            "mspace creation failed \n");

    assert(heap_mspace != 0);
    INFO(NVSHMEM_INIT, "[%d] mspace ptr: %p", state->mype, heap_mspace);

    mspace_track_large_chunks(heap_mspace, 1);

    state->heap_mspace = heap_mspace;

out:
    return status;
}

int nvshmemi_cleanup_memory_space(nvshmemi_state_t *state)
{
    destroy_mspace(state->heap_mspace);
    return 0;
}

int nvshmemi_cleanup_symmetric_heap(nvshmemi_state_t *state)
{
    INFO(NVSHMEM_INIT, "In nvshmemi_cleanup_symmetric_heap()");
    int status = 0;

    if (!state->peer_heap_base)
        goto out;

    // TODO: work required in destroying mspace
    // status = nvshmemi_cleanup_memory_space (state);
    // NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out, "memory space cleanup failed
    // \n");

    for (int i = 0; i < state->npes; i++)
    {
        if ((i == state->mype) && (state->heap_base != NULL))
        {
            for (size_t j = 0; j < state->handles.size(); j++)
            {
                status = cleanup_local_handles(
                    &state->handles[j][i * state->num_initialized_transports], state);
                NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                                      "cleanup local handles failed \n");
            }
#if CUDA_VERSION >= 11000
            if (nvshmemi_use_cuda_vmm)
            {
                for (uint32_t i = 0; i < state->cumem_handles.size(); i++)
                {
                    status = CUPFN(nvshmemi_cuda_syms, cuMemRelease(state->cumem_handles[i]));
                    NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                                          "cuMemRelease failed \n");
                }
                state->cumem_handles.clear();
            }
            else
#endif
            {
                status = cudaFree(state->peer_heap_base[i]);
                NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                                      "cudaFree failed \n");
            }

            continue;
        }

        if (state->peer_heap_base[i])
        {
            int j;
            for (j = 0; j < state->num_initialized_transports; j++)
            {
                if ((((state->transport_bitmap) & (1 << j)) &&
                     (state->transports[j]->cap[i] & NVSHMEM_TRANSPORT_CAP_MAP)) == 0)
                    continue;

                status = state->transports[j]->host_ops.unmap((state->peer_heap_base[i]),
                                                              state->heap_size);
                NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out, "unmap failed \n");

                for (size_t k = 0; k < state->handles.size(); k++)
                {
                    close(*(int *)&state->handles[k][i * state->num_initialized_transports]);
                }
            }
        }
    }
    state->handles.clear();
    state->idx_in_handles.clear();
    nvshmemi_cleanup_memory_space(state);

    free(state->peer_heap_base);
    INFO(NVSHMEM_INIT, "Leaving nvshmemi_cleanup_symmetric_heap()");

out:
    return status;
}

int nvshmemi_setup_local_heap(nvshmemi_state_t *state)
{
    int status;
    size_t alignbytes = MALLOC_ALIGNMENT;
    size_t heapextra;
    size_t tmp;

#if CUDA_VERSION >= 11000
    CUmemAllocationProp prop = {};
    prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop.location.id = static_cast<int>(state->device_id);
    prop.requestedHandleTypes = CU_MEM_HANDLE_TYPE_POSIX_FILE_DESCRIPTOR;
    prop.allocFlags.gpuDirectRDMACapable = 1;

    status = CUPFN(nvshmemi_cuda_syms,
                   cuMemGetAllocationGranularity(&cumem_granularity, &prop,
                                                 CU_MEM_ALLOC_GRANULARITY_RECOMMENDED));
    NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                          "cuMemGetAllocationGranularity failed \n");
    if (nvshmemi_use_cuda_vmm)
    {
        cumem_granularity = std::max(nvshmemi_options.CUMEM_GRANULARITY, cumem_granularity);
        cumem_granularity = cumem_granularity < NVSHMEMI_MAX_HANDLE_LENGTH
                                ? cumem_granularity
                                : NVSHMEMI_MAX_HANDLE_LENGTH;
    }
    else
    {
#endif
        cumem_granularity = NVSHMEMI_MAX_HANDLE_LENGTH;
#if CUDA_VERSION >= 11000
    }
#endif

    assert((cumem_granularity & (cumem_granularity - 1)) == 0);
    tmp = cumem_granularity;
    log2_cumem_granularity = 0;
    while (tmp >> 1)
    {
        tmp >>= 1;
        log2_cumem_granularity++;
    }

    heapextra = NUM_G_BUF_ELEMENTS * sizeof(g_elem_t) + nvshmemi_get_teams_mem_requirement() +
                G_COALESCING_BUF_SIZE + 4 * alignbytes +
                20 * alignbytes; // alignbytes, providing capacity for 2 allocations for
                                 // the library and 10 allocations for the user

    INFO(NVSHMEM_INIT, "nvshmemi_setup_local_heap, heapextra = %lld", heapextra);
    if (nvshmemi_use_cuda_vmm)
    {
        state->heap_size = std::max(nvshmemi_options.MAX_MEMORY_PER_GPU, heapextra);
    }
    else
    {
        state->heap_size = nvshmemi_options.SYMMETRIC_SIZE + heapextra;
    }
    state->heap_size =
        ((state->heap_size + cumem_granularity - 1) / cumem_granularity) * cumem_granularity;

    state->physical_heap_size = 0;
#if CUDA_VERSION >= 11000
    if (nvshmemi_use_cuda_vmm)
    {
        status = CUPFN(nvshmemi_cuda_syms,
                       cuMemAddressReserve((CUdeviceptr *)&state->global_heap_base,
                                           nvshmemi_options.MAX_P2P_GPUS * state->heap_size,
                                           alignbytes, (CUdeviceptr)NULL, 0));
        NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                              "cuMemAddressReserve failed \n");

        state->heap_base = (void *)((uintptr_t)state->global_heap_base);

        status = nvshmemi_setup_memory_space(state, state->heap_base, state->physical_heap_size);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                              "memory space initialization failed \n");
    }
    else
#endif
    {
        bool data =
            true; /*A boolean attribute which when set, ensures that synchronous memory operations
                     initiated on the region of memory that ptr points to will always synchronize.*/

        status = cudaMalloc(&state->heap_base, state->heap_size);
        NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_OUT_OF_MEMORY, out,
                              "cuMemAlloc failed \n");

        status =
            CUPFN(nvshmemi_cuda_syms, cuPointerSetAttribute(&data, CU_POINTER_ATTRIBUTE_SYNC_MEMOPS,
                                                            (CUdeviceptr)state->heap_base));
        NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_OUT_OF_MEMORY, out,
                              "cuPointerSetAttribute failed \n");

        INFO(NVSHMEM_INIT, "[%d] heap base: %p NVSHMEM_SYMMETRIC_SIZE %lu total %lu heapextra %lu",
             state->mype, state->heap_base, nvshmemi_options.SYMMETRIC_SIZE, state->heap_size,
             heapextra);

        status = nvshmemi_setup_memory_space(state, state->heap_base, state->heap_size);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                              "memory space initialization failed \n");
    }

out:
    if (status)
    {
        if (state->heap_base && nvshmemi_use_cuda_vmm == 0)
            cudaFree(state->heap_base);
    }
    return status;
}

#ifdef NVSHMEM_IBGDA_SUPPORT
int nvshmemi_gather_mem_handles(nvshmem_mem_handle_t *local_handles, nvshmemi_state_t *state,
                                uint64_t heap_offset, size_t size)
#else
int nvshmemi_gather_mem_handles(nvshmem_mem_handle_t *local_handles, nvshmemi_state_t *state)
#endif
{
    int status = nvshmemi_boot_handle.allgather(
        (void *)local_handles, (void *)(state->handles.back().data()),
        sizeof(nvshmem_mem_handle_t) * state->num_initialized_transports, &nvshmemi_boot_handle);
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                          "allgather of mem handles failed \n");

#ifdef NVSHMEM_IBGDA_SUPPORT
    for (int i = 0; i < state->num_initialized_transports; ++i)
    {
        nvshmem_transport_t tcurr = state->transports[i];
        if ((state->transport_bitmap & (1 << i)) && tcurr->host_ops.add_device_remote_mem_handles)
        {
            status = tcurr->host_ops.add_device_remote_mem_handles(
                tcurr, state->num_initialized_transports, state->handles.back().data(), heap_offset,
                size);
            NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                                  "add_device_remote_mem_handles failed \n");

            nvshmemi_gic_device_state_t *gic_state;
            nvshmemx_gic_get_device_state((void **)&gic_state);
            status = nvshmemi_gic_set_device_state(gic_state);
            NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                                  "Failed to update application ibgda state \n");
            status = nvshmemi_update_device_state();
            NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                                  "Failed to update internal ibgda state \n");
        }
    }
#endif

out:
    return status;
}

static int register_heap_handle(nvshmemi_state_t *state, nvshmem_mem_handle_t *mem_handle_in,
                                void *buf, size_t size)
{
    int status = 0;
    nvshmem_transport_t *transports = (nvshmem_transport_t *)state->transports;
    nvshmem_mem_handle_t local_handles[state->num_initialized_transports];
    nvshmem_mem_handle_t *map_handles;
    char *buf_ptr = (char *)buf;
    void *buf_map;
    // assuming symmetry of transports across all PEs
    memset(local_handles, 0, sizeof(nvshmem_mem_handle_t) * state->num_initialized_transports);

    for (int i = 0; i < state->num_initialized_transports; i++)
    {
        if ((state->transport_bitmap & (1 << i)) && transports[i]->host_ops.get_mem_handle)
        {
            INFO(NVSHMEM_INIT, "calling get_mem_handle for transport: %d buf: %p size: %lu", i, buf,
                 size);
            /* This is an unfortunate remnant of the fact that we need to map all memory at once if
             * not using VMM. */
            if (!nvshmemi_use_cuda_vmm &&
                (transports[i]->cap[state->mype] & NVSHMEM_TRANSPORT_CAP_MAP) &&
                !nvshmemi_no_vmm_heap_initialized)
            {
                status = transports[i]->host_ops.get_mem_handle(
                    (nvshmem_mem_handle_t *)(local_handles + i), mem_handle_in, state->heap_base,
                    state->heap_size, transports[i], false);
            }
            else if (nvshmemi_use_cuda_vmm ||
                     (!(transports[i]->cap[state->mype] & NVSHMEM_TRANSPORT_CAP_MAP)))
            {
                status = transports[i]->host_ops.get_mem_handle(
                    (nvshmem_mem_handle_t *)(local_handles + i), mem_handle_in, buf, size,
                    transports[i], false);
                /* This is an unfortunate remnant of the fact that we need to map all memory at once
                 * if not using VMM. */
            }
            else
            {
                map_handles = (nvshmem_mem_handle_t *)&state->handles.front();
                local_handles[i] = map_handles[i];
            }

            NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                                  "transport get memhandle failed \n");
            INFO(NVSHMEM_INIT, "[%d] get_mem_handle transport %d handles %p", state->mype, i,
                 local_handles + i);
        }
    }

    state->handles.push_back(
        vector<nvshmem_mem_handle_t>(state->num_initialized_transports * state->npes));
#ifdef NVSHMEM_IBGDA_SUPPORT
    status = nvshmemi_gather_mem_handles(local_handles, state, state->physical_heap_size, size);
#else
    status = nvshmemi_gather_mem_handles(local_handles, state);
#endif
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                          "allgather of mem handles failed \n");

#if CUDA_VERSION >= 11000
    if (nvshmemi_use_cuda_vmm)
    {
        ipcHandle *myIpcHandle = NULL;
        pid_t pid = getpid();
        IPC_CHECK(ipcOpenSocket(myIpcHandle));

        /* Wait for all processes to open their sockets */
        status = nvshmemi_boot_handle.barrier(&nvshmemi_boot_handle);

        for (std::map<pid_t, int>::iterator it1 = p2p_processes.begin(); it1 != p2p_processes.end();
             ++it1)
        {
            pid_t sending_process = it1->first;
            if (pid == sending_process)
            {
                for (std::map<pid_t, int>::iterator it2 = p2p_processes.begin();
                     it2 != p2p_processes.end(); ++it2)
                {
                    pid_t receiving_process = it2->first;
                    if (pid != receiving_process)
                    { /* Dont sent to yourself */
                        IPC_CHECK(ipcSendFd(myIpcHandle, *(int *)local_handles, receiving_process));
                    }
                }
            }
            else
            {
                IPC_CHECK(ipcRecvFd(myIpcHandle,
                                    (int *)&state->handles
                                        .back()[it1->second * state->num_initialized_transports]));
            }
            status = nvshmemi_boot_handle.barrier(&nvshmemi_boot_handle);
        }
        IPC_CHECK(ipcCloseSocket(myIpcHandle));
        mspace_add_new_chunk(state->heap_mspace,
                             (char *)state->heap_base + state->physical_heap_size, size);
    }
#endif
    for (int i = ((state->mype + 1) % state->npes); i != state->mype; i = ((i + 1) % state->npes))
    {
        for (int j = 0; j < state->num_initialized_transports; j++)
        {
            if ((state->transport_map[state->mype * state->npes + i] & (1 << j)) &&
                (transports[j]->cap[i] & NVSHMEM_TRANSPORT_CAP_MAP))
            {
                if (!nvshmemi_use_cuda_vmm && !nvshmemi_no_vmm_heap_initialized)
                {
                    status = transports[j]->host_ops.map(
                        (state->peer_heap_base + i), state->heap_size,
                        &state->handles.back()[i * state->num_initialized_transports + j]);
                }
                else if (nvshmemi_use_cuda_vmm)
                {
                    buf_map = (void *)(buf_ptr - (char *)state->heap_base +
                                       (char *)(state->peer_heap_base[i]));
                    status = transports[j]->host_ops.map(
                        &buf_map, size,
                        &state->handles.back()[i * state->num_initialized_transports + j]);
                }
                else
                {
                    status = 0;
                }
                if (status)
                {
                    // map operation failed, remove cap of transport
                    state->transports[j]->cap[i] ^= NVSHMEM_TRANSPORT_CAP_MAP;
                    status = 0;
                    continue;
                }
                char *hex = nvshmemu_hexdump(
                    &state->handles.back()[i * state->num_initialized_transports + j],
                    sizeof(CUipcMemHandle));
                INFO(NVSHMEM_INIT, "[%d] cuIpcOpenMemHandle fromhandle 0x%s", state->mype, hex);
                free(hex);

                INFO(NVSHMEM_INIT, "[%d] cuIpcOpenMemHandle tobuf %p", state->mype,
                     *(state->peer_heap_base + i));
                break;
            }
        }
    }

    /*
     * Don't set nvshmemi_no_vmm_heap_initialized until
     * we have completed the above for loop and
     * properly mapped the heap from each PE.
     */
    if (!nvshmemi_use_cuda_vmm)
    {
        nvshmemi_no_vmm_heap_initialized = true;
    }

    for (size_t i = 0; i < (size / cumem_granularity); i++)
    {
        nvshmemi_state->idx_in_handles.push_back(
            make_tuple(nvshmemi_state->handles.size() - 1,
                       (char *)state->heap_base + state->physical_heap_size, size));
    }
    state->physical_heap_size += size;

out:
    return status;
}

static int register_heap_handles(nvshmemi_state_t *state, nvshmem_mem_handle_t *mem_handle_in,
                                 void *buf, size_t size)
{
    int status = 0;
    size_t remaining_size;
    size_t registration_size;
    char *buf_calculation = (char *)buf;

    if (size == 0)
    {
        return NVSHMEMX_ERROR_INVALID_VALUE;
    }

    remaining_size = size;
    do
    {
        registration_size = size;

        status =
            register_heap_handle(state, mem_handle_in, (void *)buf_calculation, registration_size);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                              "register heap handle failed \n");

        assert(remaining_size >= registration_size);
        remaining_size -= registration_size;
        buf_calculation += registration_size;
    } while (remaining_size);

out:
    return status;
}

int nvshmemi_setup_symmetric_heap(nvshmemi_state_t *state)
{
    int status;
    int p2p_counter = 1;
    pid_t pid = 0;
    uint64_t myHostHash;
    pid_t *peer_pids = NULL;
    uint64_t *hostHash = NULL;
    nvshmem_transport_t *transports = (nvshmem_transport_t *)state->transports;

    nvshmemi_no_vmm_heap_initialized = false;
    state->peer_heap_base_actual = (void **)calloc(state->npes, sizeof(void *));
    NVSHMEMI_NULL_ERROR_JMP(state->peer_heap_base_actual, status, NVSHMEMX_ERROR_OUT_OF_MEMORY, out,
                            "failed allocating space for peer heap base \n");

    status = nvshmemi_boot_handle.allgather((void *)&state->heap_base,
                                            (void *)state->peer_heap_base_actual, sizeof(void *),
                                            &nvshmemi_boot_handle);
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                          "allgather of heap base ptrsmem handle failed \n");

    state->peer_heap_base = (void **)calloc(state->npes, sizeof(void *));
    NVSHMEMI_NULL_ERROR_JMP(state->peer_heap_base, status, NVSHMEMX_ERROR_OUT_OF_MEMORY, out,
                            "failed allocating space for peer heap base \n");

    state->peer_heap_base[state->mype] = state->heap_base;
    if (!nvshmemi_use_cuda_vmm)
    {
        status = register_heap_handles(state, NULL, state->heap_base, state->heap_size);
        NVSHMEMI_NE_ERROR_JMP(status, 0, NVSHMEMX_ERROR_INTERNAL, out,
                              "register heap handle failed \n");
    }
    else
    {
        for (int i = ((state->mype + 1) % state->npes); i != state->mype;
             i = ((i + 1) % state->npes))
        {
            for (int j = 0; j < state->num_initialized_transports; j++)
            {
                if (state->transport_map[state->mype * state->npes + i] & (1 << j))
                {
                    if (transports[j]->cap[i] & NVSHMEM_TRANSPORT_CAP_MAP)
                    {
                        state->peer_heap_base[i] = (void *)((uintptr_t)state->global_heap_base +
                                                            state->heap_size * p2p_counter++);
                        break;
                    }
                }
            }
        }

        /* Build p2p_processes that is used during dynaminic heap management */
        pid = getpid();
        peer_pids = (pid_t *)malloc(sizeof(pid_t) * state->npes);
        status = nvshmemi_boot_handle.allgather((void *)&pid, (void *)peer_pids, sizeof(pid_t),
                                                &nvshmemi_boot_handle);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out, "allgather of pids failed \n");

        myHostHash = getHostHash();
        hostHash = (uint64_t *)malloc(sizeof(uint64_t) * state->npes);
        status = nvshmemi_boot_handle.allgather((void *)&myHostHash, (void *)hostHash,
                                                sizeof(uint64_t), &nvshmemi_boot_handle);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                              "allgather of host hashes failed \n");

        for (int pe = 0; pe < state->npes; pe++)
        {
            if (hostHash[pe] == myHostHash)
            {
                p2p_processes[peer_pids[pe]] = pe;
            }
        }
        INFO(NVSHMEM_MEM, "I am connected to %lu p2p processes (including myself)",
             p2p_processes.size());
    }

out:
    if (peer_pids)
    {
        free(peer_pids);
    }
    if (hostHash)
    {
        free(hostHash);
    }
    if (status)
    {
        // if handles has been allocated, try and cleanup all heap state
        // else cleanup local handles only
        nvshmemi_cleanup_symmetric_heap(state);
        if (state->heap_base)
            cudaFree(state->heap_base);
    }
    return status;
}

#if CUDA_VERSION >= 11000
/* This function is an exclusive agent of the function nvshmemi_add_physical_memory
 * As such, there are two assumptions inherent to this function -
 * the size argument passed by nvshmemi_add_physical_memory must be aligned to cumem_granularity.
 * this function will not change the value of size.
 */
static int add_physical_memory_aligned(nvshmemi_state_t *state, size_t size)
{
    CUmemGenericAllocationHandle cumem_handle;
    CUmemAllocationProp prop = {};
    CUmemAccessDesc access;
    char *buf_calc;
    size_t remaining_size;
    size_t register_size;
    size_t adjusted_max_handle_len;
    int status;

    prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop.location.id = static_cast<int>(state->device_id);
    prop.requestedHandleTypes = CU_MEM_HANDLE_TYPE_POSIX_FILE_DESCRIPTOR;
    prop.allocFlags.gpuDirectRDMACapable = 1;

    access.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    access.location.id = state->device_id;
    access.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;

    assert(size % cumem_granularity == 0);
    assert(cumem_granularity <= NVSHMEMI_MAX_HANDLE_LENGTH);

    /* Round Down */
    adjusted_max_handle_len = cumem_granularity * (NVSHMEMI_MAX_HANDLE_LENGTH / cumem_granularity);

    remaining_size = size;
    do
    {
        buf_calc = (char *)state->heap_base + state->physical_heap_size;
        register_size =
            remaining_size > adjusted_max_handle_len ? adjusted_max_handle_len : remaining_size;

        status = CUPFN(nvshmemi_cuda_syms, cuMemCreate(&cumem_handle, register_size,
                                                       (const CUmemAllocationProp *)&prop, 0));
        NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                              "cuMemCreate failed \n");
        state->cumem_handles.push_back(cumem_handle);

        status = CUPFN(nvshmemi_cuda_syms,
                       cuMemMap((CUdeviceptr)buf_calc, register_size, 0, cumem_handle, 0));
        NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                              "cuMemMap failed \n");

        status = CUPFN(nvshmemi_cuda_syms, cuMemSetAccess((CUdeviceptr)buf_calc, register_size,
                                                          (const CUmemAccessDesc *)&access, 1));
        NVSHMEMI_NE_ERROR_JMP(status, CUDA_SUCCESS, NVSHMEMX_ERROR_INTERNAL, out,
                              "cuMemSetAccess failed \n");

        status = register_heap_handle(state, (nvshmem_mem_handle_t *)&cumem_handle, buf_calc,
                                      register_size);
        NVSHMEMI_NE_ERROR_JMP(status, 0, NVSHMEMX_ERROR_INTERNAL, out,
                              "register heap handle failed \n");

        remaining_size -= register_size;
    } while (remaining_size > 0);
    /* The above loop works iff remaining_size and adjusted_max_handle_len are both multiples of
     * cumem_granularity */
out:
    return status;
}

int nvshmemi_add_physical_memory(size_t size)
{
    int status;
    nvshmemi_state_t *state = nvshmemi_state;

    size = ((size + cumem_granularity - 1) / cumem_granularity) * cumem_granularity;
    INFO(NVSHMEM_MEM, "Adding new physical backing of size %zu bytes", size);

    status = add_physical_memory_aligned(state, size);
    NVSHMEMI_NE_ERROR_JMP(status, 0, NVSHMEMX_ERROR_INTERNAL, out, "add physical memory failed \n");

    status = nvshmemi_boot_handle.barrier(
        &nvshmemi_boot_handle); /* Wait for all PEs to setup the new memory */
out:
    if (status)
    {
        nvshmemi_cleanup_symmetric_heap(state);
        if (state->heap_base)
            cudaFree(state->heap_base);
    }
    return status;
}
#endif

extern "C"
{
    void *nvshmemi_malloc(size_t size)
    {
        int status = 0;
        void *ptr = NULL;

        status = check_for_symmetry(size);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INVALID_VALUE, out,
                              "symmetry check for size failed\n");

        ptr = mspace_malloc(nvshmemi_state->heap_mspace, size);
#if CUDA_VERSION >= 11000
        if (nvshmemi_use_cuda_vmm)
        {
            if ((size > 0) && (ptr == NULL))
            {
                nvshmemi_add_physical_memory(size);
                ptr = mspace_malloc(nvshmemi_state->heap_mspace, size);
            }
            goto update_device_state;
        }
        else
#endif
        {
            if ((size > 0) && (ptr == NULL))
            {
                NVSHMEMI_ERROR_EXIT(
                    "nvshmem malloc failed (hint: check if total allocation has exceeded NVSHMEM "
                    "symmetric size = %zu, NVSHMEM symmetric size can be increased using "
                    "NVSHMEM_SYMMETRIC_SIZE environment variable) \n",
                    nvshmemi_options.SYMMETRIC_SIZE);
            }
        }

    update_device_state:
        status = nvshmemi_update_device_state();
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                              "nvshmemi_update_device_state failed\n");

        INFO(NVSHMEM_INIT, "[%d] allocated %zu bytes from mspace: %p ptr: %p", nvshmemi_state->mype,
             size, nvshmemi_state->heap_mspace, ptr);

    out:
        if (status)
        {
            if (ptr)
            {
                mspace_free(nvshmemi_state->heap_mspace, ptr);
                ptr = NULL;
            }
        }
        return ptr;
    }
}

void *nvshmemi_calloc(size_t count, size_t size)
{
    int status = 0;
    void *ptr = NULL;

    status = check_for_symmetry(size);
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INVALID_VALUE, out,
                          "symmetry check for size failed\n");

    ptr = mspace_calloc(nvshmemi_state->heap_mspace, count, size);
#if CUDA_VERSION >= 11000
    if (nvshmemi_use_cuda_vmm)
    {
        if ((size > 0) && (ptr == NULL))
        {
            nvshmemi_add_physical_memory(size);
            ptr = mspace_calloc(nvshmemi_state->heap_mspace, count, size);
        }
        goto update_device_state;
    }
    else
#endif
    {
        if (size > 0 && count > 0 && ptr == NULL)
        {
            NVSHMEMI_ERROR_EXIT("nvshmem calloc failed \n");
        }
    }

update_device_state:
    status = nvshmemi_update_device_state();
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                          "nvshmemi_update_device_state failed\n");

    INFO(NVSHMEM_INIT, "[%d] calloc allocated %zu bytes from mspace: %p ptr: %p",
         nvshmemi_state->mype, size, nvshmemi_state->heap_mspace, ptr);
out:
    if (status)
    {
        if (ptr)
        {
            mspace_free(nvshmemi_state->heap_mspace, ptr);
            ptr = NULL;
        }
    }
    return ptr;
}

void *nvshmem_malloc(size_t size)
{
    void *ptr = NULL;

    NVTX_FUNC_RANGE_IN_GROUP(ALLOC);

    NVSHMEMU_THREAD_CS_ENTER();
    (*nvshmemi_check_state_and_init_fn_ptr)();

    ptr = nvshmemi_malloc(size);

    nvshmemi_barrier_all();

    NVSHMEMU_THREAD_CS_EXIT();

    return ptr;
}

void *nvshmem_calloc(size_t count, size_t size)
{
    void *ptr = NULL;

    NVTX_FUNC_RANGE_IN_GROUP(ALLOC);

    NVSHMEMU_THREAD_CS_ENTER();
    (*nvshmemi_check_state_and_init_fn_ptr)();

    ptr = nvshmemi_calloc(count, size);

    nvshmemi_barrier_all();

    NVSHMEMU_THREAD_CS_EXIT();

    return ptr;
}

void *nvshmemi_align(size_t alignment, size_t size)
{
    void *ptr = NULL;
    int status = 0;

    status = check_for_symmetry(size);
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INVALID_VALUE, out,
                          "symmetry check for size failed\n");

    ptr = mspace_memalign(nvshmemi_state->heap_mspace, alignment, size);
#if CUDA_VERSION >= 11000
    if (nvshmemi_use_cuda_vmm)
    {
        if ((size > 0) && (ptr == NULL))
        {
            nvshmemi_add_physical_memory(size + alignment);
            ptr = mspace_memalign(nvshmemi_state->heap_mspace, alignment, size);
        }
        goto update_device_state;
    }
    else
#endif
    {
        if ((size > 0) && (ptr == NULL))
        {
            NVSHMEMI_ERROR_EXIT("nvshmem align failed \n");
        }
    }

update_device_state:
    status = nvshmemi_update_device_state();
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out,
                          "nvshmemi_update_device_state failed\n");

out:
    if (status)
    {
        if (ptr)
        {
            mspace_free(nvshmemi_state->heap_mspace, ptr);
            ptr = NULL;
        }
    }
    return ptr;
}

void *nvshmem_align(size_t alignment, size_t size)
{
    void *ptr = NULL;

    NVTX_FUNC_RANGE_IN_GROUP(ALLOC);

    NVSHMEMU_THREAD_CS_ENTER();
    (*nvshmemi_check_state_and_init_fn_ptr)();

    ptr = nvshmemi_align(alignment, size);

    nvshmemi_barrier_all();

    NVSHMEMU_THREAD_CS_EXIT();

    return ptr;
}

void nvshmemi_free(void *ptr)
{
    if (ptr == NULL)
        return;

    INFO(NVSHMEM_INIT, "[%d] freeing buf: %p", nvshmemi_state->mype, ptr);

    mspace_free(nvshmemi_state->heap_mspace, ptr);

    nvshmemi_update_device_state();
}

void nvshmem_free(void *ptr)
{
    NVTX_FUNC_RANGE_IN_GROUP(ALLOC);

    NVSHMEMU_THREAD_CS_ENTER();

    NVSHMEMI_CHECK_INIT_STATUS();

    nvshmemi_barrier_all();

    nvshmemi_free(ptr);

    NVSHMEMU_THREAD_CS_EXIT();
}

void *nvshmem_ptr(const void *ptr, int pe)
{
    if (ptr >= nvshmemi_state->heap_base)
    {
        uintptr_t offset = (char *)ptr - (char *)nvshmemi_state->heap_base;

        if (offset < nvshmemi_state->heap_size)
        {
            void *peer_addr = nvshmemi_state->peer_heap_base[pe];
            if (peer_addr != NULL)
                peer_addr = (void *)((char *)peer_addr + offset);
            return peer_addr;
        }
    }

    return NULL;
}

/* for now, pick up the first transport - in the future, we need to only support one remote
 * transport. */
static struct nvshmem_transport *nvshmemi_get_remote_transport()
{
    struct nvshmem_transport *t = NULL;

    for (int j = 0; j < nvshmemi_state->num_initialized_transports; j++)
    {
        if (nvshmemi_state->transport_bitmap & (1 << j))
        {
            t = nvshmemi_state->transports[j];
        }
    }

    return t;
}

static int buffer_register(nvshmem_transport_t transport, void *addr, size_t length)
{
    nvshmem_local_buf_cache_t *cache = (nvshmem_local_buf_cache_t *)transport->cache_handle;
    nvshmem_local_buf_handle_t *handle;
    size_t i;
    size_t selected_index;
    size_t number_of_handles;
    size_t remaining_length;
    void *heap_end;
    int status = 0;
    int lock_status = EBUSY;
    cudaPointerAttributes attr;
#if CUDART_VERSION < 11000
    bool register_with_cuda = false;
#endif
#ifdef NVSHMEM_IBGDA_SUPPORT
    nvshmemi_gic_device_state_t *gic_state;
#endif

    if (length == 0)
    {
        NVSHMEMI_ERROR_PRINT("Unable to register zero length buffers.\n");
        return NVSHMEMX_ERROR_INVALID_VALUE;
    }

    number_of_handles = (length + NVSHMEMI_MAX_HANDLE_LENGTH - 1) / NVSHMEMI_MAX_HANDLE_LENGTH;
    status = cudaPointerGetAttributes(&attr, addr);
#if CUDART_VERSION >= 11000
    if (status != cudaSuccess)
    {
        NVSHMEMI_ERROR_PRINT("Unable to query pointer attributes.\n");
        /* clear CUDA error string. */
        cudaGetLastError();
        return NVSHMEMX_ERROR_INTERNAL;
    }
#else
    if (status != cudaSuccess)
    {
        /* clear CUDA error string. */
        cudaGetLastError();
        if (status == cudaErrorInvalidValue)
        {
            register_with_cuda = true;
        }
        else
        {
            NVSHMEMI_ERROR_PRINT("Unable to query pointer attributes.\n");
            return NVSHMEMX_ERROR_INTERNAL;
        }
    }
#endif

    if (attr.type == cudaMemoryTypeManaged)
    {
        NVSHMEMI_ERROR_PRINT("Unable to register managed memory as it can migrate.\n");
        return NVSHMEMX_ERROR_INVALID_VALUE;
    }

    heap_end = (void *)((char *)nvshmemi_state->heap_base + nvshmemi_state->heap_size);
    if (addr >= nvshmemi_state->heap_base && addr < heap_end)
    {
        NVSHMEMI_ERROR_PRINT(
            "Unable to register nvshmem heap memory. It is registered by default.\n");
        return NVSHMEMX_ERROR_INVALID_VALUE;
    }

    handle =
        (nvshmem_local_buf_handle_t *)calloc(number_of_handles, sizeof(nvshmem_local_buf_handle_t));
    if (handle == NULL)
    {
        NVSHMEMI_ERROR_PRINT("Unable to resize the registered buffer array.\n");
        return NVSHMEMX_ERROR_OUT_OF_MEMORY;
    }

    if (transport)
    {
        for (i = 0; i < number_of_handles; i++)
        {
            handle[i].handle = (nvshmem_mem_handle_t *)calloc(1, sizeof(nvshmem_mem_handle_t));
            handle[i].linked_with_prev = true;
            handle[i].linked_with_next = true;
            if (handle[i].handle == NULL)
            {
                NVSHMEMI_ERROR_PRINT("Unable to resize the registered buffer array.\n");
                status = NVSHMEMX_ERROR_OUT_OF_MEMORY;
                goto out_error_unlocked;
            }
        }
        handle[0].linked_with_prev = false;
        handle[number_of_handles - 1].linked_with_next = false;
    }

    while (lock_status == EBUSY)
    {
        lock_status = pthread_rwlock_wrlock(&cache->buffer_lock);
    }

    if (lock_status != 0)
    {
        NVSHMEMI_ERROR_PRINT("Unable to acquire buffer registration lock with errno %d\n",
                             lock_status);
        status = NVSHMEMX_ERROR_INTERNAL;
        goto out_error_unlocked;
    }

#if CUDART_VERSION >= 11000
    if (attr.type == cudaMemoryTypeUnregistered)
    {
#else
    if (register_with_cuda)
    {
#endif
        if (!nvshmemi_state->host_memory_registration_supported)
        {
            NVSHMEMI_ERROR_PRINT(
                "Unable to register host memory for this device as it doesn't support UVA.\n");
            status = NVSHMEMX_ERROR_INVALID_VALUE;
            goto out_unlock;
        }
        status = cudaHostRegister(addr, length, cudaHostRegisterDefault);
        if (status)
        {
            NVSHMEMI_ERROR_PRINT("Unable to register host memory with CUDA.\n");
            status = NVSHMEMX_ERROR_INTERNAL;
            goto out_unlock;
        }
        for (i = 0; i < number_of_handles; i++)
        {
            handle[i].registered_by_us = true;
        }
    }

    /* We only need to register unregistered host buffers if there is no remote transport.
     * CUDA memory and registered host memory are already mapped into the address space
     * so nothing to register.
     */
    if (transport == NULL && handle->registered_by_us == false)
    {
        status = 0;
        free(handle);
        goto out_unlock;
    }

    if ((cache->array_used + number_of_handles) >= cache->array_size)
    {
        size_t new_array_size;
        void *new_buf;

        if (number_of_handles > cache->array_size)
        {
            new_array_size = cache->array_size + number_of_handles;
        }
        else
        {
            new_array_size = cache->array_size * 2;
        }

        assert(new_array_size < (SIZE_MAX / sizeof(nvshmem_local_buf_handle_t)));
        new_buf = realloc(cache->buffers, new_array_size * sizeof(nvshmem_local_buf_handle_t *));
        if (new_buf == NULL)
        {
            NVSHMEMI_ERROR_PRINT("Unable to resize the registered buffer array.\n");
            status = NVSHMEMX_ERROR_OUT_OF_MEMORY;
            goto out_unlock;
        }
        cache->buffers = (nvshmem_local_buf_handle_t **)new_buf;
        cache->array_size = new_array_size;
    }

    for (i = 0; i < cache->array_used; i++)
    {
        nvshmem_local_buf_handle_t *tmp_handle = cache->buffers[i];
        if (addr > tmp_handle->ptr)
        {
            void *max_addr;
            max_addr = (void *)((char *)tmp_handle->ptr + tmp_handle->length);
            if (addr < max_addr)
            {
                NVSHMEMI_ERROR_PRINT("Unable to register overlapping memory regions.\n");
                status = NVSHMEMX_ERROR_INVALID_VALUE;
                goto out_unlock;
            }
            continue;
        }
        else if (addr == tmp_handle->ptr)
        {
            NVSHMEMI_ERROR_PRINT("Unable to double register memory regions.\n");
            status = NVSHMEMX_ERROR_INVALID_VALUE;
            goto out_unlock;
            /* addr < tmp_handle->ptr */
        }
        else
        {
            break;
        }
    }

    selected_index = i;
    remaining_length = length;
    for (i = 0; i < number_of_handles; i++)
    {
        char *addr_calc = (char *)addr;
        size_t register_length = remaining_length > NVSHMEMI_MAX_HANDLE_LENGTH
                                     ? NVSHMEMI_MAX_HANDLE_LENGTH
                                     : remaining_length;
        if (transport)
        {
            assert(register_length < NVSHMEMI_DMA_BUF_MAX_LENGTH);
            status = transport->host_ops.get_mem_handle(handle[i].handle, NULL, (void *)addr_calc,
                                                        register_length, transport, true);
            if (status)
            {
                NVSHMEMI_ERROR_PRINT("Unable to assign new memory handle.\n");
                goto out_unlock;
            }
        }
        handle[i].ptr = (void *)addr_calc;
        handle[i].length = register_length;

        addr_calc += register_length;
        remaining_length -= register_length;
    }

    assert(selected_index < cache->array_size);
    if (selected_index < cache->array_used)
    {
        memmove(&cache->buffers[selected_index + number_of_handles],
                &cache->buffers[selected_index],
                sizeof(nvshmem_local_buf_handle_t *) * (cache->array_used - selected_index));
    }
    for (i = 0; i < number_of_handles; i++)
    {
        cache->buffers[selected_index + i] = &handle[i];
        cache->array_used++;
    }
#ifdef NVSHMEM_IBGDA_SUPPORT
    if (transport->type == NVSHMEM_TRANSPORT_LIB_CODE_IBGDA)
    {
        nvshmemx_gic_get_device_state((void **)&gic_state);
        status = nvshmemi_gic_set_device_state(gic_state);
        NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out_unlock,
                              "Failed to update application ibgda state \n");
    }
#endif
    status = nvshmemi_update_device_state();
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out_unlock,
                          "nvshmemi_update_device_state failed\n");

out_unlock:
    pthread_rwlock_unlock(&cache->buffer_lock);
    if (status == 0)
    {
        return 0;
    }

out_error_unlocked:
    if (handle->registered_by_us)
    {
        cudaHostUnregister(addr);
    }
    if (handle)
    {
        for (i = 0; i < number_of_handles; i++)
        {
            if (handle[i].handle)
            {
                free(handle[i].handle);
            }
        }
        free(handle);
    }
    return status;
}

int nvshmemx_buffer_register(void *addr, size_t length)
{
    int status_global = NVSHMEMX_SUCCESS;
    int status_local;

    for (int i = 0; i < nvshmemi_state->num_initialized_transports; i++)
    {
        if (nvshmemi_state->transport_bitmap & (1 << i))
        {
            status_local = buffer_register(nvshmemi_state->transports[i], addr, length);
            if (status_local)
            {
                NVSHMEMI_ERROR_PRINT("Buffer addition for transport %d failed with error %d\n", i,
                                     status_local);
                status_global = status_local;
            }
        }
    }

    return status_global;
}

static int buffer_unregister(nvshmem_transport_t transport, void *addr)
{
    nvshmem_local_buf_cache_t *cache = (nvshmem_local_buf_cache_t *)transport->cache_handle;
    size_t i;
    size_t num_mem_handles = 0;
    int lock_status = EBUSY;
    bool linked_with_next;
    int status = 0;
#ifdef NVSHMEM_IBGDA_SUPPORT
    nvshmemi_gic_device_state_t *gic_state;
#endif

    while (lock_status == EBUSY)
    {
        lock_status = pthread_rwlock_wrlock(&cache->buffer_lock);
    }

    if (lock_status != 0)
    {
        NVSHMEMI_ERROR_PRINT("Unable to acquire buffer registration lock with errno %d\n",
                             lock_status);
        return NVSHMEMX_ERROR_INTERNAL;
    }

    for (i = 0; i < cache->array_used; i++)
    {
        nvshmem_local_buf_handle_t *tmp_handle = cache->buffers[i];
        nvshmem_local_buf_handle_t *tmp_handle_next;
        if (addr > tmp_handle->ptr)
        {
            continue;
        }
        else if (addr == tmp_handle->ptr)
        {
            do
            {
                linked_with_next = tmp_handle->linked_with_next;
                if (transport)
                {
                    transport->host_ops.release_mem_handle(tmp_handle->handle, transport);
                    free(tmp_handle->handle);
                }

                if (tmp_handle->registered_by_us && !tmp_handle->linked_with_prev)
                {
                    cudaHostUnregister(tmp_handle->ptr);
                }
                num_mem_handles++;
                tmp_handle_next = tmp_handle + 1;
                tmp_handle = tmp_handle_next;
            } while (linked_with_next);
            free(cache->buffers[i]);

            if ((i + num_mem_handles) < cache->array_used)
            {
                memmove(&cache->buffers[i], &cache->buffers[i + num_mem_handles],
                        sizeof(nvshmem_local_buf_handle_t *) *
                            (cache->array_used - (i + num_mem_handles - 1)));
            }

            cache->array_used -= num_mem_handles;
#ifdef NVSHMEM_IBGDA_SUPPORT
            if (transport->type == NVSHMEM_TRANSPORT_LIB_CODE_IBGDA)
            {
                nvshmemx_gic_get_device_state((void **)&gic_state);
                status = nvshmemi_gic_set_device_state(gic_state);
                NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out_unlock,
                                      "Failed to update application ibgda state \n");
            }
#endif
            status = nvshmemi_update_device_state();
            NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMX_ERROR_INTERNAL, out_unlock,
                                  "nvshmemi_update_device_state failed\n");
            goto out_unlock;
            /* addr < tmp_handle->ptr*/
        }
        else
        {
            break;
        }
    }

    status = NVSHMEMX_ERROR_INVALID_VALUE;
out_unlock:
    pthread_rwlock_unlock(&cache->buffer_lock);
    return status;
}

int nvshmemx_buffer_unregister(void *addr)
{
    int status_global = NVSHMEMX_SUCCESS;
    int status_local;

    for (int i = 0; i < nvshmemi_state->num_initialized_transports; i++)
    {
        if (nvshmemi_state->transport_bitmap & (1 << i))
        {
            status_local = buffer_unregister(nvshmemi_state->transports[i], addr);
            if (status_local)
            {
                NVSHMEMI_ERROR_PRINT("Buffer removal for transport %d failed with error %d\n", i,
                                     status_local);
                status_global = status_local;
            }
        }
    }

    return status_global;
}

static void buffer_unregister_all(nvshmem_transport_t transport)
{
    nvshmem_local_buf_cache_t *cache = (nvshmem_local_buf_cache_t *)transport->cache_handle;
    int lock_status = EBUSY;

    while (lock_status == EBUSY)
    {
        lock_status = pthread_rwlock_wrlock(&cache->buffer_lock);
    }

    if (lock_status != 0)
    {
        NVSHMEMI_ERROR_PRINT(
            "Unable to acquire buffer registration lock with errno %d. Unregister all function "
            "failed.\n",
            lock_status);
        return;
    }

    for (size_t i = 0; i < cache->array_used; i++)
    {
        if (transport)
        {
            transport->host_ops.release_mem_handle(cache->buffers[i]->handle, transport);
        }

        if (cache->buffers[i]->registered_by_us && !cache->buffers[i]->linked_with_prev)
        {
            cudaHostUnregister(cache->buffers[i]->ptr);
        }

        free(cache->buffers[i]->handle);
        free(cache->buffers[i]);
    }

    cache->array_used = 0;
    pthread_rwlock_unlock(&cache->buffer_lock);

    return;
}

void nvshmemx_buffer_unregister_all()
{
    for (int i = 0; i < nvshmemi_state->num_initialized_transports; i++)
    {
        if (nvshmemi_state->transport_bitmap & (1 << i))
        {
            buffer_unregister_all(nvshmemi_state->transports[i]);
        }
    }

    return;
}

struct nvshmem_mem_handle *nvshmemi_get_registered_buffer_handle(nvshmem_transport_t transport,
                                                                 void *addr, size_t *len)
{
    nvshmem_local_buf_cache_t *cache = (nvshmem_local_buf_cache_t *)transport->cache_handle;
    nvshmem_local_buf_handle_t *tmp_handle;
    size_t min, max, mid;
    void *max_addr;
    size_t max_len;
    int lock_status = EBUSY;
    struct nvshmem_mem_handle *ret_handle = NULL;

    while (lock_status == EBUSY)
    {
        lock_status = pthread_rwlock_rdlock(&cache->buffer_lock);
    }

    if (lock_status != 0)
    {
        NVSHMEMI_ERROR_PRINT("Unable to acquire buffer registration lock with errno %d.\n",
                             lock_status);
        return ret_handle;
    }

    if (cache->array_used == 0)
    {
        goto out_unlock;
    }

    min = 0;
    max = cache->array_used;
    do
    {
        mid = (max - min) / 2 + min;
        /* We have gone past the end of the loop. */
        if (mid >= cache->array_used)
        {
            break;
        }
        tmp_handle = cache->buffers[mid];
        if (addr > tmp_handle->ptr)
        {
            max_addr = (void *)((char *)tmp_handle->ptr + tmp_handle->length);
            max_len = (uint64_t)((char *)max_addr - (char *)addr);
            if (addr < max_addr)
            {
                *len = *len < max_len ? *len : max_len;
                ret_handle = tmp_handle->handle;
                goto out_unlock;
            }
            min = mid + 1;
        }
        else if (addr == tmp_handle->ptr)
        {
            *len = *len < tmp_handle->length ? *len : tmp_handle->length;
            ret_handle = tmp_handle->handle;
            goto out_unlock;
        }
        else
        {
            if (mid == 0)
            {
                break;
            }
            max = mid - 1;
        }
    } while (max >= min);

    NVSHMEMI_ERROR_PRINT("Unable to find a reference to the requested buffer address.\n");

out_unlock:
    pthread_rwlock_unlock(&cache->buffer_lock);
    return ret_handle;
}

void nvshmemi_local_mem_cache_fini(nvshmem_transport_t transport,
                                   nvshmem_local_buf_cache_t *cache)
{
    buffer_unregister_all(transport);
    pthread_rwlock_destroy(&cache->buffer_lock);
    free(cache->buffers);
    free(cache);
}

int nvshmemi_local_mem_cache_init(nvshmem_local_buf_cache_t **cache)
{
    nvshmem_local_buf_cache_t *cache_ptr;
    int status;

    cache_ptr = (nvshmem_local_buf_cache_t *)calloc(1, sizeof(nvshmem_local_buf_cache_t));
    NVSHMEMI_NULL_ERROR_JMP(cache_ptr, status, NVSHMEMI_INTERNAL_ERROR, err,
                            "Unable to allocate cache.\n");
    cache_ptr->buffers = (nvshmem_local_buf_handle_t **)calloc(
        NVSHMEMI_LOCAL_BUF_CACHE_DEFAULT_SIZE, sizeof(nvshmem_local_buf_handle_t *));
    cache_ptr->array_size = NVSHMEMI_LOCAL_BUF_CACHE_DEFAULT_SIZE;
    cache_ptr->array_used = 0;
    status = pthread_rwlock_init(&cache_ptr->buffer_lock, NULL);
    NVSHMEMI_NZ_ERROR_JMP(status, NVSHMEMI_INTERNAL_ERROR, err, "Unable to init cache lock.\n");

    *cache = cache_ptr;
    return NVSHMEMI_SUCCESS;
err:
    free(cache_ptr->buffers);
    free(cache_ptr);
    return status;
}
