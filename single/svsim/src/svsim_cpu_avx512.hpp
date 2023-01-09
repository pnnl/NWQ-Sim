// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: svsim_cpu_avx512.hpp
// AVX512 implementation for SV-Sim CPU backends.
// ---------------------------------------------------------------------------
//
#ifndef SVSIM_CPU_AVX512_HPP_
#define SVSIM_CPU_AVX512_HPP_
#include <immintrin.h>

/***********************************************
 * Key Macros
 ***********************************************/
//Common
#define PUT(arr,i,val) (_mm512_i64scatter_pd((arr),(i),(val),8))
#define GET(arr,i) (_mm512_i64gather_pd((i),(arr),8))
#define BARR while(0){}
//Load 1st qubit coefficient
#define LOAD_Q0 \
    const __m512d el0_real=GET(sv_real,pos0);\
    const __m512d el0_imag=GET(sv_imag,pos0);\
    const __m512d el1_real=GET(sv_real,pos1);\
    const __m512d el1_imag=GET(sv_imag,pos1);
//Save 1st qubit coefficient
#define STORE_Q0 \
    PUT(sv_real,pos0,sv_real_pos0);\
    PUT(sv_imag,pos0,sv_imag_pos0);\
    PUT(sv_real,pos1,sv_real_pos1);\
    PUT(sv_imag,pos1,sv_imag_pos1);
//Load 2nd qubit coefficient
#define LOAD_Q1 \
    const __m512d el2_real=GET(sv_real,pos2);\
    const __m512d el2_imag=GET(sv_imag,pos2);\
    const __m512d el3_real=GET(sv_real,pos3);\
    const __m512d el3_imag=GET(sv_imag,pos3);
//Save 2nd qubit coefficient
#define STORE_Q1 \
    PUT(sv_real,pos2,sv_real_pos2);\
    PUT(sv_imag,pos2,sv_imag_pos2);\
    PUT(sv_real,pos3,sv_real_pos3);\
    PUT(sv_imag,pos3,sv_imag_pos3);
//macro for C2/C4 gates
#define DIV2E(x,y) _mm512_srli_epi64((x),(y))
#define MOD2E(x,y) _mm512_and_epi64((x),_mm512_set1_epi64(((IdxType)1<<(y))-(IdxType)1))
#define EXP2E(x) _mm512_set1_epi64((IdxType)1<<(x))
#define EXP2S(x) ((IdxType)1<<(x))
#define SV4IDX(x) _mm512_set1_epi64( ((x>>1)&1)*EXP2S(ctrl) + ((x&1)*EXP2S(qubit)) )
#define SV16IDX(x) _mm512_set1_epi64((((x)>>3)&1)*EXP2S(qubit0)+(((x)>>2)&1)*EXP2S(qubit1)+(((x)>>1)&1)*EXP2S(qubit2)+(((x)&1)*EXP2S(qubit3)))
//Define MG-BSP machine operation header (Optimized version)
#define OP_HEAD \
    const __m512i cons0 = _mm512_set1_epi64((sim->half_dim)-1); \
    const __m512i cons1 = _mm512_set1_epi64(((IdxType)1<<qubit)-1); \
    const __m512i cons2 = _mm512_set1_epi64((IdxType)1<<qubit); \
    const __m512i base =_mm512_set_epi64(0,1,2,3,4,5,6,7); \
    _Pragma("omp for schedule(auto)") \
    for (IdxType i=0; i<(sim->half_dim);i+=8){ \
        __m512i idx=_mm512_set1_epi64(i); \
        idx=_mm512_add_epi64(idx,base); \
        __m512i tmp = _mm512_and_epi64(idx,cons0); \
        __m512i outer = _mm512_srli_epi64(tmp,qubit);\
        __m512i inner = _mm512_and_epi64(idx,cons1); \
        __m512i offset = _mm512_slli_epi64(outer,qubit+1); \
        __m512i pos0 = _mm512_add_epi64(offset,inner); \
        __m512i pos1 = _mm512_add_epi64(pos0,cons2); 
//Define operation header for 2-qubit
#define OP_HEAD_2Q \
    assert(ctrl != qubit); \
    const IdxType p = min(ctrl, qubit); \
    const IdxType q = max(ctrl, qubit); \
    const __m512i base = _mm512_set_epi64(0,1,2,3,4,5,6,7); \
    _Pragma("omp for schedule(auto)") \
    for (IdxType i=0; i<((sim->dim)>>2); i+=8) \
    { \
        __m512i idx=_mm512_add_epi64(_mm512_set1_epi64(i),base); \
        const __m512i term0 = MOD2E(idx,p); \
        const __m512i term1 = _mm512_mullo_epi64(MOD2E(DIV2E(idx,p),q-p-1),EXP2E(p+1)); \
        const __m512i term2 = _mm512_mullo_epi64(DIV2E(DIV2E(idx,p),q-p-1),EXP2E(q+1)); \
        const __m512i term = _mm512_add_epi64(term2,_mm512_add_epi64(term1,term0)); \
        __m512i pos0  = _mm512_add_epi64(term,SV4IDX(0)); \
        __m512i pos1  = _mm512_add_epi64(term,SV4IDX(1)); \
        __m512i pos2  = _mm512_add_epi64(term,SV4IDX(2)); \
        __m512i pos3  = _mm512_add_epi64(term,SV4IDX(3));

//Define operation tail
#define OP_TAIL } BARR;

/***********************************************
 * Gate Implementation
 ***********************************************/

//============== C1 Gate ================
void C1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 0]), el0_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 0]), el0_imag)), 
                                         _mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 1]), el1_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 1]), el1_imag)));
    __m512d sv_imag_pos0 = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 0]), el0_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 0]), el0_real)), 
                                         _mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 1]), el1_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 1]), el1_real)));
    __m512d sv_real_pos1 = _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 2]), el0_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 2]), el0_imag)), 
                                         _mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 3]), el1_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 3]), el1_imag)));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 2]), el0_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 2]), el0_real)), 
                                         _mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 3]), el1_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 3]), el1_real)));
    STORE_Q0;
    OP_TAIL;
}


//============== C2 Gate ================
//Arbitrary 2-qubit gate
inline void C2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag,
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const __m512d el0_real = GET(sv_real, pos0);
    const __m512d el0_imag = GET(sv_imag, pos0);
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    const __m512d el2_real = GET(sv_real, pos2);
    const __m512d el2_imag = GET(sv_imag, pos2);
    const __m512d el3_real = GET(sv_real, pos3);
    const __m512d el3_imag = GET(sv_imag, pos3);
    //Real part
    __m512d sv_real_pos0 = _mm512_add_pd( _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 0]), el0_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 0]), el0_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 1]), el1_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 1]), el1_imag)))),
                                          _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 2]), el2_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 2]), el2_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 3]), el3_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 3]), el3_imag)))));
    __m512d sv_real_pos1 = _mm512_add_pd( _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 4]), el0_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 4]), el0_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 5]), el1_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 5]), el1_imag)))),
                                          _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 6]), el2_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 6]), el2_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 7]), el3_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 7]), el3_imag)))));
    __m512d sv_real_pos2 = _mm512_add_pd( _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 8]), el0_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 8]), el0_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 9]), el1_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 9]), el1_imag)))),
                                          _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[10]), el2_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[10]), el2_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[11]), el3_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[11]), el3_imag)))));
    __m512d sv_real_pos3 = _mm512_add_pd( _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[12]), el0_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[12]), el0_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[13]), el1_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[13]), el1_imag)))),
                                          _mm512_add_pd(_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[14]), el2_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[14]), el2_imag)), 
                                                       (_mm512_sub_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[15]), el3_real), _mm512_mul_pd(_mm512_set1_pd(gm_imag[15]), el3_imag)))));
    //Imag part
    __m512d sv_imag_pos0 = _mm512_add_pd( _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 0]), el0_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 0]), el0_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 1]), el1_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 1]), el1_real)))),
                                          _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 2]), el2_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 2]), el2_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 3]), el3_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 3]), el3_real)))));
    __m512d sv_imag_pos1 = _mm512_add_pd( _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 4]), el0_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 4]), el0_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 5]), el1_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 5]), el1_real)))),
                                          _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 6]), el2_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 6]), el2_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 7]), el3_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 7]), el3_real)))));
    __m512d sv_imag_pos2 = _mm512_add_pd( _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 8]), el0_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 8]), el0_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[ 9]), el1_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[ 9]), el1_real)))),
                                          _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[10]), el2_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[10]), el2_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[11]), el3_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[11]), el3_real)))));
    __m512d sv_imag_pos3 = _mm512_add_pd( _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[12]), el0_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[12]), el0_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[13]), el1_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[13]), el1_real)))),
                                          _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[14]), el2_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[14]), el2_real)), 
                                                       (_mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(gm_real[15]), el3_imag), _mm512_mul_pd(_mm512_set1_pd(gm_imag[15]), el3_real)))));
    PUT(sv_real, pos0, sv_real_pos0);
    PUT(sv_imag, pos0, sv_imag_pos0);
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    PUT(sv_real, pos2, sv_real_pos2);
    PUT(sv_imag, pos2, sv_imag_pos2);
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//#define SV16IDX(x) _mm512_add_epi64(_mm512_add_epi64( _mm512_mullo_epi64(_mm512_and_epi64(_mm512_srli_epi64((x),3),1),EXP2E(qubit0)),\
//_mm512_mullo_epi64(_mm512_and_epi64(_mm512_srli_epi64((x),2),1),EXP2E(qubit1))),\
//_mm512_add_epi64( _mm512_mullo_epi64(_mm512_and_epi64(_mm512_srli_epi64((x),1),1),EXP2E(qubit2)),\
//_mm512_mullo_epi64(_mm512_and_epi64(_mm512_srli_epi64((x),0),1),EXP2E(qubit3))))

//============== C4 Gate ================
//Arbitrary 4-qubit gate
inline void C4_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag,
        const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    assert (qubit0 != qubit1); //Non-cloning
    assert (qubit0 != qubit2); //Non-cloning
    assert (qubit0 != qubit3); //Non-cloning
    assert (qubit1 != qubit2); //Non-cloning
    assert (qubit1 != qubit3); //Non-cloning
    assert (qubit2 != qubit3); //Non-cloning
    //need to sort qubits: min->max: p, q, r, s
    const IdxType v0 = min(qubit0, qubit1);
    const IdxType v1 = min(qubit2, qubit3);
    const IdxType v2 = max(qubit0, qubit1);
    const IdxType v3 = max(qubit2, qubit3);
    const IdxType p = min(v0,v1); 
    const IdxType q = min(min(v2,v3),max(v0,v1)); 
    const IdxType r = max(min(v2,v3),max(v0,v1)); 
    const IdxType s = max(v2,v3);
    const __m512i base = _mm512_set_epi64(0,1,2,3,4,5,6,7);
    _Pragma("omp for schedule(auto)") 
    for (IdxType i=0; i<((sim->dim)>>4); i+=8)
    {
        __m512i idx=_mm512_add_epi64(_mm512_set1_epi64(i),base);
        const __m512i term0 = MOD2E(idx,p);
        const __m512i term1 = _mm512_mullo_epi64(MOD2E(DIV2E(idx,p),q-p-1),EXP2E(p+1));
        const __m512i term2 = _mm512_mullo_epi64(MOD2E(DIV2E(DIV2E(idx,p),q-p-1),r-q-1),EXP2E(q+1));
        const __m512i term3 = _mm512_mullo_epi64(MOD2E(DIV2E(DIV2E(DIV2E(idx,p),q-p-1),r-q-1),s-r-1),EXP2E(r+1));
        const __m512i term4 = _mm512_mullo_epi64(DIV2E(DIV2E(DIV2E(DIV2E(idx,p),q-p-1),r-q-1),s-r-1),EXP2E(s+1));
        const __m512i term = _mm512_add_epi64(term4,_mm512_add_epi64(term3,_mm512_add_epi64(term2,_mm512_add_epi64(term1,term0))));

        const __m512d el0_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(0)));
        const __m512d el1_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(1)));
        const __m512d el2_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(2)));
        const __m512d el3_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(3)));
        const __m512d el4_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(4)));
        const __m512d el5_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(5)));
        const __m512d el6_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(6)));
        const __m512d el7_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(7)));
        const __m512d el8_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(8)));
        const __m512d el9_real  = GET(sv_real,_mm512_add_epi64(term,SV16IDX(9)));
        const __m512d el10_real = GET(sv_real,_mm512_add_epi64(term,SV16IDX(10)));
        const __m512d el11_real = GET(sv_real,_mm512_add_epi64(term,SV16IDX(11)));
        const __m512d el12_real = GET(sv_real,_mm512_add_epi64(term,SV16IDX(12)));
        const __m512d el13_real = GET(sv_real,_mm512_add_epi64(term,SV16IDX(13)));
        const __m512d el14_real = GET(sv_real,_mm512_add_epi64(term,SV16IDX(14)));
        const __m512d el15_real = GET(sv_real,_mm512_add_epi64(term,SV16IDX(15)));

        const __m512d el0_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(0)));
        const __m512d el1_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(1)));
        const __m512d el2_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(2)));
        const __m512d el3_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(3)));
        const __m512d el4_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(4)));
        const __m512d el5_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(5)));
        const __m512d el6_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(6)));
        const __m512d el7_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(7)));
        const __m512d el8_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(8)));
        const __m512d el9_imag  = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(9)));
        const __m512d el10_imag = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(10)));
        const __m512d el11_imag = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(11)));
        const __m512d el12_imag = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(12)));
        const __m512d el13_imag = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(13)));
        const __m512d el14_imag = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(14)));
        const __m512d el15_imag = GET(sv_imag,_mm512_add_epi64(term,SV16IDX(15)));

        for (unsigned j=0; j<16; j++)
        {
            __m512d res_real = _mm512_set1_pd(0);
            __m512d res_imag = _mm512_set1_pd(0);
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el0_real,_mm512_set1_pd(gm_real[j*16+0])),
                                                             _mm512_mul_pd(el0_imag,_mm512_set1_pd(gm_imag[j*16+0]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el0_real,_mm512_set1_pd(gm_imag[j*16+0])),
                                                             _mm512_mul_pd(el0_imag,_mm512_set1_pd(gm_real[j*16+0]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el1_real,_mm512_set1_pd(gm_real[j*16+1])),
                                                             _mm512_mul_pd(el1_imag,_mm512_set1_pd(gm_imag[j*16+1]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el1_real,_mm512_set1_pd(gm_imag[j*16+1])),
                                                             _mm512_mul_pd(el1_imag,_mm512_set1_pd(gm_real[j*16+1]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el2_real,_mm512_set1_pd(gm_real[j*16+2])),
                                                             _mm512_mul_pd(el2_imag,_mm512_set1_pd(gm_imag[j*16+2]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el2_real,_mm512_set1_pd(gm_imag[j*16+2])),
                                                             _mm512_mul_pd(el2_imag,_mm512_set1_pd(gm_real[j*16+2]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el3_real,_mm512_set1_pd(gm_real[j*16+3])),
                                                             _mm512_mul_pd(el3_imag,_mm512_set1_pd(gm_imag[j*16+3]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el3_real,_mm512_set1_pd(gm_imag[j*16+3])),
                                                             _mm512_mul_pd(el3_imag,_mm512_set1_pd(gm_real[j*16+3]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el4_real,_mm512_set1_pd(gm_real[j*16+4])),
                                                             _mm512_mul_pd(el4_imag,_mm512_set1_pd(gm_imag[j*16+4]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el4_real,_mm512_set1_pd(gm_imag[j*16+4])),
                                                             _mm512_mul_pd(el4_imag,_mm512_set1_pd(gm_real[j*16+4]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el5_real,_mm512_set1_pd(gm_real[j*16+5])),
                                                             _mm512_mul_pd(el5_imag,_mm512_set1_pd(gm_imag[j*16+5]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el5_real,_mm512_set1_pd(gm_imag[j*16+5])),
                                                             _mm512_mul_pd(el5_imag,_mm512_set1_pd(gm_real[j*16+5]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el6_real,_mm512_set1_pd(gm_real[j*16+6])),
                                                             _mm512_mul_pd(el6_imag,_mm512_set1_pd(gm_imag[j*16+6]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el6_real,_mm512_set1_pd(gm_imag[j*16+6])),
                                                             _mm512_mul_pd(el6_imag,_mm512_set1_pd(gm_real[j*16+6]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el7_real,_mm512_set1_pd(gm_real[j*16+7])),
                                                             _mm512_mul_pd(el7_imag,_mm512_set1_pd(gm_imag[j*16+7]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el7_real,_mm512_set1_pd(gm_imag[j*16+7])),
                                                             _mm512_mul_pd(el7_imag,_mm512_set1_pd(gm_real[j*16+7]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el8_real,_mm512_set1_pd(gm_real[j*16+8])),
                                                             _mm512_mul_pd(el8_imag,_mm512_set1_pd(gm_imag[j*16+8]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el8_real,_mm512_set1_pd(gm_imag[j*16+8])),
                                                             _mm512_mul_pd(el8_imag,_mm512_set1_pd(gm_real[j*16+8]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el9_real,_mm512_set1_pd(gm_real[j*16+9])),
                                                             _mm512_mul_pd(el9_imag,_mm512_set1_pd(gm_imag[j*16+9]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el9_real,_mm512_set1_pd(gm_imag[j*16+9])),
                                                             _mm512_mul_pd(el9_imag,_mm512_set1_pd(gm_real[j*16+9]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el10_real,_mm512_set1_pd(gm_real[j*16+10])),
                                                             _mm512_mul_pd(el10_imag,_mm512_set1_pd(gm_imag[j*16+10]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el10_real,_mm512_set1_pd(gm_imag[j*16+10])),
                                                             _mm512_mul_pd(el10_imag,_mm512_set1_pd(gm_real[j*16+10]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el11_real,_mm512_set1_pd(gm_real[j*16+11])),
                                                             _mm512_mul_pd(el11_imag,_mm512_set1_pd(gm_imag[j*16+11]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el11_real,_mm512_set1_pd(gm_imag[j*16+11])),
                                                             _mm512_mul_pd(el11_imag,_mm512_set1_pd(gm_real[j*16+11]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el12_real,_mm512_set1_pd(gm_real[j*16+12])),
                                                             _mm512_mul_pd(el12_imag,_mm512_set1_pd(gm_imag[j*16+12]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el12_real,_mm512_set1_pd(gm_imag[j*16+12])),
                                                             _mm512_mul_pd(el12_imag,_mm512_set1_pd(gm_real[j*16+12]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el13_real,_mm512_set1_pd(gm_real[j*16+13])),
                                                             _mm512_mul_pd(el13_imag,_mm512_set1_pd(gm_imag[j*16+13]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el13_real,_mm512_set1_pd(gm_imag[j*16+13])),
                                                             _mm512_mul_pd(el13_imag,_mm512_set1_pd(gm_real[j*16+13]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el14_real,_mm512_set1_pd(gm_real[j*16+14])),
                                                             _mm512_mul_pd(el14_imag,_mm512_set1_pd(gm_imag[j*16+14]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el14_real,_mm512_set1_pd(gm_imag[j*16+14])),
                                                             _mm512_mul_pd(el14_imag,_mm512_set1_pd(gm_real[j*16+14]))));
            res_real = _mm512_add_pd(res_real, _mm512_sub_pd(_mm512_mul_pd(el15_real,_mm512_set1_pd(gm_real[j*16+15])),
                                                             _mm512_mul_pd(el15_imag,_mm512_set1_pd(gm_imag[j*16+15]))));
            res_imag = _mm512_add_pd(res_imag, _mm512_add_pd(_mm512_mul_pd(el15_real,_mm512_set1_pd(gm_imag[j*16+15])),
                                                             _mm512_mul_pd(el15_imag,_mm512_set1_pd(gm_real[j*16+15]))));
            PUT(sv_real, _mm512_add_epi64(term, SV16IDX(j)), res_real);
            PUT(sv_imag, _mm512_add_epi64(term, SV16IDX(j)), res_imag);
        }
    }
    BARR;
}
#endif
