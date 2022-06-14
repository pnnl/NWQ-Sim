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
// File: dmsim_cpu_avx512.hpp
// AVX512 implementation for DM-Sim CPU backends.
// ---------------------------------------------------------------------------
//
#ifndef NWQSIM_CPU_AVX512_HPP_
#define NWQSIM_CPU_AVX512_HPP_
#include <immintrin.h>

/***********************************************
 * Key Macros
 ***********************************************/
//Common
#define PUT(arr,i,val) (_mm512_i64scatter_pd((arr),(i),(val),8))
#define GET(arr,i) (_mm512_i64gather_pd((i),(arr),8))
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
#define OP_TAIL } _Pragma("omp barrier")  

/***********************************************
 * Gate Implementation
 ***********************************************/
//============== X Gate ================
inline void X_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = el1_real;
    __m512d sv_imag_pos0 = el1_imag;
    __m512d sv_real_pos1 = el0_real;
    __m512d sv_imag_pos1 = el0_imag;
    STORE_Q0;
    OP_TAIL;
}
//============== Y Gate ================
inline void Y_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 =  el1_imag; 
    __m512d sv_imag_pos0 = -el1_real;
    __m512d sv_real_pos1 = -el0_imag;
    __m512d sv_imag_pos1 =  el0_real;
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateY Gate ================
inline void ConjugateY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = -el1_imag; 
    __m512d sv_imag_pos0 =  el1_real;
    __m512d sv_real_pos1 =  el0_imag;
    __m512d sv_imag_pos1 = -el0_real;
    STORE_Q0;
    OP_TAIL;
}
//============== Z Gate ================
inline void Z_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 = -el1_real;
    __m512d sv_imag_pos1 = -el1_imag;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== H Gate ================
inline void H_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    const __m512d s2i_v = _mm512_set1_pd(S2I);
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el0_real, el1_real)); 
    __m512d sv_imag_pos0 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el0_imag, el1_imag)); 
    __m512d sv_real_pos1 = _mm512_mul_pd(s2i_v, _mm512_sub_pd(el0_real, el1_real)); 
    __m512d sv_imag_pos1 = _mm512_mul_pd(s2i_v, _mm512_sub_pd(el0_imag, el1_imag)); 
    STORE_Q0;
    OP_TAIL;
}
//============== S Gate ================
inline void S_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 = -el1_imag;
    __m512d sv_imag_pos1 =  el1_real;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== SDG Gate ================
inline void SDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    __m512d el1_real = GET(sv_real, pos1);
    __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 =  el1_imag;
    __m512d sv_imag_pos1 = -el1_real;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== T Gate ================
inline void T_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    const __m512d s2i_v = _mm512_set1_pd(S2I);
    OP_HEAD;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 = _mm512_mul_pd(s2i_v, _mm512_sub_pd(el1_real, el1_imag)); 
    __m512d sv_imag_pos1 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el1_real, el1_imag)); 
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== TDG Gate ================
inline void TDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    const __m512d s2i_v = _mm512_set1_pd(S2I);
    OP_HEAD;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el1_real, el1_imag)); 
    __m512d sv_imag_pos1 = _mm512_mul_pd(s2i_v, _mm512_add_pd(-el1_real, el1_imag)); 
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== RI Gate ================
inline void RI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d ri_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d ri_imag = _mm512_set1_pd(-sin(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd(_mm512_mul_pd(el0_real, ri_real), _mm512_mul_pd(el0_imag, ri_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd(_mm512_mul_pd(el0_real, ri_imag), _mm512_mul_pd(el0_imag, ri_real));
    __m512d sv_real_pos1 = _mm512_sub_pd(_mm512_mul_pd(el1_real, ri_real), _mm512_mul_pd(el1_imag, ri_imag));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(el1_real, ri_imag), _mm512_mul_pd(el1_imag, ri_real));
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateRI Gate ================
inline void ConjugateRI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d ri_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d ri_imag = _mm512_set1_pd(sin(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd(_mm512_mul_pd(el0_real, ri_real), _mm512_mul_pd(el0_imag, ri_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd(_mm512_mul_pd(el0_real, ri_imag), _mm512_mul_pd(el0_imag, ri_real));
    __m512d sv_real_pos1 = _mm512_sub_pd(_mm512_mul_pd(el1_real, ri_real), _mm512_mul_pd(el1_imag, ri_imag));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(el1_real, ri_imag), _mm512_mul_pd(el1_imag, ri_real));
    STORE_Q0;
    OP_TAIL;
}
//============== RX Gate ================
inline void RX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d rx_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d rx_imag = _mm512_set1_pd(-sin(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd( _mm512_mul_pd(rx_real, el0_real), _mm512_mul_pd(rx_imag, el1_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd( _mm512_mul_pd(rx_real, el0_imag), _mm512_mul_pd(rx_imag, el1_real));
    __m512d sv_real_pos1 = _mm512_add_pd(-_mm512_mul_pd(rx_imag, el0_imag), _mm512_mul_pd(rx_real, el1_real));
    __m512d sv_imag_pos1 = _mm512_add_pd( _mm512_mul_pd(rx_imag, el0_real), _mm512_mul_pd(rx_real, el1_imag));
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateRX Gate ================
inline void ConjugateRX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d rx_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d rx_imag = _mm512_set1_pd(sin(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd( _mm512_mul_pd(rx_real, el0_real), _mm512_mul_pd(rx_imag, el1_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd( _mm512_mul_pd(rx_real, el0_imag), _mm512_mul_pd(rx_imag, el1_real));
    __m512d sv_real_pos1 = _mm512_add_pd(-_mm512_mul_pd(rx_imag, el0_imag), _mm512_mul_pd(rx_real, el1_real));
    __m512d sv_imag_pos1 = _mm512_add_pd( _mm512_mul_pd(rx_imag, el0_real), _mm512_mul_pd(rx_real, el1_imag));
    STORE_Q0;
    OP_TAIL;
}
//============== RY Gate ================
inline void RY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e1_real = _mm512_set1_pd(-sin(HALF*theta));
    const __m512d e2_real = _mm512_set1_pd(sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_add_pd(_mm512_mul_pd(e0_real, el0_real), _mm512_mul_pd(e1_real, el1_real));
    __m512d sv_imag_pos0 = _mm512_add_pd(_mm512_mul_pd(e0_real, el0_imag), _mm512_mul_pd(e1_real, el1_imag));
    __m512d sv_real_pos1 = _mm512_add_pd(_mm512_mul_pd(e2_real, el0_real), _mm512_mul_pd(e3_real, el1_real));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(e2_real, el0_imag), _mm512_mul_pd(e3_real, el1_imag));
    STORE_Q0;
    OP_TAIL;
}
//============== RZ Gate ================
inline void RZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e0_imag = _mm512_set1_pd(-sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(sin(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd(_mm512_mul_pd(el0_real, e0_real), _mm512_mul_pd(el0_imag, e0_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd(_mm512_mul_pd(el0_real, e0_imag), _mm512_mul_pd(el0_imag, e0_real));
    __m512d sv_real_pos1 = _mm512_sub_pd(_mm512_mul_pd(el1_real, e3_real), _mm512_mul_pd(el1_imag, e3_imag));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(el1_real, e3_imag), _mm512_mul_pd(el1_imag, e3_real));
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateRZ Gate ================
inline void ConjugateRZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e0_imag = _mm512_set1_pd(sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(-sin(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd(_mm512_mul_pd(el0_real, e0_real), _mm512_mul_pd(el0_imag, e0_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd(_mm512_mul_pd(el0_real, e0_imag), _mm512_mul_pd(el0_imag, e0_real));
    __m512d sv_real_pos1 = _mm512_sub_pd(_mm512_mul_pd(el1_real, e3_real), _mm512_mul_pd(el1_imag, e3_imag));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(el1_real, e3_imag), _mm512_mul_pd(el1_imag, e3_real));
    STORE_Q0;
    OP_TAIL;
}
//============== SX Gate ================
inline void SX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    const __m512d half_v = _mm512_set1_pd(HALF);
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el0_real,el0_imag),_mm512_add_pd(el1_real,el1_imag)));
    __m512d sv_imag_pos0 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el0_imag,el0_real),_mm512_sub_pd(el1_imag,el1_real)));
    __m512d sv_real_pos1 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el0_real,el0_imag),_mm512_sub_pd(el1_real,el1_imag)));
    __m512d sv_imag_pos1 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el0_imag,el0_real),_mm512_add_pd(el1_imag,el1_real)));
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateSX Gate ================
inline void ConjugateSX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    const __m512d half_v = _mm512_set1_pd(HALF);
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el0_real,el0_imag),_mm512_sub_pd(el1_real,el1_imag)));
    __m512d sv_imag_pos0 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el0_imag,el0_real),_mm512_add_pd(el1_imag,el1_real)));
    __m512d sv_real_pos1 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el0_real,el0_imag),_mm512_add_pd(el1_real,el1_imag)));
    __m512d sv_imag_pos1 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el0_imag,el0_real),_mm512_sub_pd(el1_imag,el1_real)));
    STORE_Q0;
    OP_TAIL;
}

//============== P Gate ================
inline void P_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    const __m512d e3_real = _mm512_set1_pd(cos(theta));
    const __m512d e3_imag = _mm512_set1_pd(sin(theta));
    OP_HEAD;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 = _mm512_sub_pd(_mm512_mul_pd(e3_real, el1_real), _mm512_mul_pd(e3_imag, el1_imag));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(e3_real, el1_imag), _mm512_mul_pd(e3_imag, el1_real));
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}

//============== ConjugateP Gate ================
inline void ConjugateP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const ValType theta, const IdxType qubit)
{
    const __m512d e3_real = _mm512_set1_pd(cos(theta));
    const __m512d e3_imag = _mm512_set1_pd(-sin(theta));
    OP_HEAD;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos1 = _mm512_sub_pd(_mm512_mul_pd(e3_real, el1_real), _mm512_mul_pd(e3_imag, el1_imag));
    __m512d sv_imag_pos1 = _mm512_add_pd(_mm512_mul_pd(e3_real, el1_imag), _mm512_mul_pd(e3_imag, el1_real));
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}

//============== U Gate ================
inline void U_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const ValType phi, const ValType lambda, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e1_real = _mm512_set1_pd(-cos(lambda)*sin(HALF*theta));
    const __m512d e1_imag = _mm512_set1_pd(-sin(lambda)*sin(HALF*theta));
    const __m512d e2_real = _mm512_set1_pd(cos(phi)*sin(HALF*theta));
    const __m512d e2_imag = _mm512_set1_pd(sin(phi)*sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(phi+lambda)*cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(sin(phi+lambda)*cos(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el0_real), _mm512_mul_pd(e1_real, el1_real)), _mm512_mul_pd(e1_imag, el1_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el0_imag), _mm512_mul_pd(e1_real, el1_imag)), _mm512_mul_pd(e1_imag, el1_real));
    __m512d sv_real_pos1 = _mm512_add_pd( _mm512_sub_pd(_mm512_mul_pd(e2_real, el0_real), _mm512_mul_pd(e2_imag, el0_imag)), 
                                          _mm512_sub_pd(_mm512_mul_pd(e3_real, el1_real), _mm512_mul_pd(e3_imag, el1_imag))); 
    __m512d sv_imag_pos1 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e2_real, el0_imag), _mm512_mul_pd(e2_imag, el0_real)), 
                                          _mm512_add_pd(_mm512_mul_pd(e3_real, el1_imag), _mm512_mul_pd(e3_imag, el1_real))); 
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateU Gate ================
inline void ConjugateU_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const ValType phi, const ValType lambda, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e1_real = _mm512_set1_pd(-cos(lambda)*sin(HALF*theta));
    const __m512d e1_imag = _mm512_set1_pd(sin(lambda)*sin(HALF*theta));
    const __m512d e2_real = _mm512_set1_pd(cos(phi)*sin(HALF*theta));
    const __m512d e2_imag = _mm512_set1_pd(-sin(phi)*sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(phi+lambda)*cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(-sin(phi+lambda)*cos(HALF*theta));
    OP_HEAD;
    LOAD_Q0;
    __m512d sv_real_pos0 = _mm512_sub_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el0_real), _mm512_mul_pd(e1_real, el1_real)), _mm512_mul_pd(e1_imag, el1_imag));
    __m512d sv_imag_pos0 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el0_imag), _mm512_mul_pd(e1_real, el1_imag)), _mm512_mul_pd(e1_imag, el1_real));
    __m512d sv_real_pos1 = _mm512_add_pd( _mm512_sub_pd(_mm512_mul_pd(e2_real, el0_real), _mm512_mul_pd(e2_imag, el0_imag)), 
                                          _mm512_sub_pd(_mm512_mul_pd(e3_real, el1_real), _mm512_mul_pd(e3_imag, el1_imag))); 
    __m512d sv_imag_pos1 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e2_real, el0_imag), _mm512_mul_pd(e2_imag, el0_real)), 
                                          _mm512_add_pd(_mm512_mul_pd(e3_real, el1_imag), _mm512_mul_pd(e3_imag, el1_real))); 
    STORE_Q0;
    OP_TAIL;
}
//============== CX Gate ================
inline void CX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = el3_real;
    __m512d sv_imag_pos2 = el3_imag;
    __m512d sv_real_pos3 = el2_real;
    __m512d sv_imag_pos3 = el2_imag;
    STORE_Q1;
    OP_TAIL;
}
//============== CY Gate ================
inline void CY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 =  el3_imag; 
    __m512d sv_imag_pos2 = -el3_real;
    __m512d sv_real_pos3 = -el2_imag;
    __m512d sv_imag_pos3 =  el2_real;
    STORE_Q1;
    OP_TAIL;
}
//============== Conjugate CY Gate ================
inline void ConjugateCY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = -el3_imag; 
    __m512d sv_imag_pos2 =  el3_real;
    __m512d sv_real_pos3 =  el2_imag;
    __m512d sv_imag_pos3 = -el2_real;
    STORE_Q1;
    OP_TAIL;
}
//============== CZ Gate ================
inline void CZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos3);
    const __m512d el3_imag = GET(sv_imag, pos3);
    __m512d sv_real_pos3 = -el3_real;
    __m512d sv_imag_pos3 = -el3_imag;
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CH Gate ================
inline void CH_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d s2i_v = _mm512_set1_pd(S2I);
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el2_real, el3_real)); 
    __m512d sv_imag_pos2 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el2_imag, el3_imag)); 
    __m512d sv_real_pos3 = _mm512_mul_pd(s2i_v, _mm512_sub_pd(el2_real, el3_real)); 
    __m512d sv_imag_pos3 = _mm512_mul_pd(s2i_v, _mm512_sub_pd(el2_imag, el3_imag)); 
    STORE_Q1;
    OP_TAIL;
}
//============== CS Gate ================
inline void CS_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos3);
    const __m512d el3_imag = GET(sv_imag, pos3);
    __m512d sv_real_pos3 = -el3_imag;
    __m512d sv_imag_pos3 =  el3_real;
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CSDG Gate ================
inline void CSDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos3);
    const __m512d el3_imag = GET(sv_imag, pos3);
    __m512d sv_real_pos3 =  el3_imag;
    __m512d sv_imag_pos3 = -el3_real;
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CT Gate ================
inline void CT_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d s2i_v = _mm512_set1_pd(S2I);
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos3);
    const __m512d el3_imag = GET(sv_imag, pos3);
    __m512d sv_real_pos3 = _mm512_mul_pd(s2i_v, _mm512_sub_pd(el3_real, el3_imag)); 
    __m512d sv_imag_pos3 = _mm512_mul_pd(s2i_v, _mm512_add_pd(el3_real, el3_imag)); 
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CTDG Gate ================
inline void CTDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d s2i_v = _mm512_set1_pd(S2I);
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos3);
    const __m512d el3_imag = GET(sv_imag, pos3);
    __m512d sv_real_pos3 = _mm512_mul_pd(s2i_v, _mm512_add_pd( el3_real, el3_imag)); 
    __m512d sv_imag_pos3 = _mm512_mul_pd(s2i_v, _mm512_add_pd(-el3_real, el3_imag)); 
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CRI Gate ================
inline void CRI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d ri_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d ri_imag = _mm512_set1_pd(-sin(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd(_mm512_mul_pd(el2_real, ri_real), _mm512_mul_pd(el2_imag, ri_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd(_mm512_mul_pd(el2_real, ri_imag), _mm512_mul_pd(el2_imag, ri_real));
    __m512d sv_real_pos3 = _mm512_sub_pd(_mm512_mul_pd(el3_real, ri_real), _mm512_mul_pd(el3_imag, ri_imag));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(el3_real, ri_imag), _mm512_mul_pd(el3_imag, ri_real));
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCRI Gate ================
inline void ConjugateCRI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d ri_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d ri_imag = _mm512_set1_pd(sin(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd(_mm512_mul_pd(el2_real, ri_real), _mm512_mul_pd(el2_imag, ri_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd(_mm512_mul_pd(el2_real, ri_imag), _mm512_mul_pd(el2_imag, ri_real));
    __m512d sv_real_pos3 = _mm512_sub_pd(_mm512_mul_pd(el3_real, ri_real), _mm512_mul_pd(el3_imag, ri_imag));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(el3_real, ri_imag), _mm512_mul_pd(el3_imag, ri_real));
    STORE_Q1;
    OP_TAIL;
}
//============== CRX Gate ================
inline void CRX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d rx_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d rx_imag = _mm512_set1_pd(-sin(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd( _mm512_mul_pd(rx_real, el2_real), _mm512_mul_pd(rx_imag, el3_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd( _mm512_mul_pd(rx_real, el2_imag), _mm512_mul_pd(rx_imag, el3_real));
    __m512d sv_real_pos3 = _mm512_add_pd(-_mm512_mul_pd(rx_imag, el2_imag), _mm512_mul_pd(rx_real, el3_real));
    __m512d sv_imag_pos3 = _mm512_add_pd( _mm512_mul_pd(rx_imag, el2_real), _mm512_mul_pd(rx_real, el3_imag));
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCRX Gate ================
inline void ConjugateCRX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d rx_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d rx_imag = _mm512_set1_pd(sin(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd( _mm512_mul_pd(rx_real, el2_real), _mm512_mul_pd(rx_imag, el3_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd( _mm512_mul_pd(rx_real, el2_imag), _mm512_mul_pd(rx_imag, el3_real));
    __m512d sv_real_pos3 = _mm512_add_pd(-_mm512_mul_pd(rx_imag, el2_imag), _mm512_mul_pd(rx_real, el3_real));
    __m512d sv_imag_pos3 = _mm512_add_pd( _mm512_mul_pd(rx_imag, el2_real), _mm512_mul_pd(rx_real, el3_imag));
    STORE_Q1;
    OP_TAIL;
}
//============== CRY Gate ================
inline void CRY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e1_real = _mm512_set1_pd(-sin(HALF*theta));
    const __m512d e2_real = _mm512_set1_pd(sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_add_pd(_mm512_mul_pd(e0_real, el2_real), _mm512_mul_pd(e1_real, el3_real));
    __m512d sv_imag_pos2 = _mm512_add_pd(_mm512_mul_pd(e0_real, el2_imag), _mm512_mul_pd(e1_real, el3_imag));
    __m512d sv_real_pos3 = _mm512_add_pd(_mm512_mul_pd(e2_real, el2_real), _mm512_mul_pd(e3_real, el3_real));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(e2_real, el2_imag), _mm512_mul_pd(e3_real, el3_imag));
    STORE_Q1;
    OP_TAIL;
}
//============== CRZ Gate ================
inline void CRZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e0_imag = _mm512_set1_pd(-sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(sin(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd(_mm512_mul_pd(el2_real, e0_real), _mm512_mul_pd(el2_imag, e0_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd(_mm512_mul_pd(el2_real, e0_imag), _mm512_mul_pd(el2_imag, e0_real));
    __m512d sv_real_pos3 = _mm512_sub_pd(_mm512_mul_pd(el3_real, e3_real), _mm512_mul_pd(el3_imag, e3_imag));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(el3_real, e3_imag), _mm512_mul_pd(el3_imag, e3_real));
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCRZ Gate ================
inline void ConjugateCRZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e0_imag = _mm512_set1_pd(sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(-sin(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd(_mm512_mul_pd(el2_real, e0_real), _mm512_mul_pd(el2_imag, e0_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd(_mm512_mul_pd(el2_real, e0_imag), _mm512_mul_pd(el2_imag, e0_real));
    __m512d sv_real_pos3 = _mm512_sub_pd(_mm512_mul_pd(el3_real, e3_real), _mm512_mul_pd(el3_imag, e3_imag));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(el3_real, e3_imag), _mm512_mul_pd(el3_imag, e3_real));
    STORE_Q1;
    OP_TAIL;
}
//============== CSX Gate ================
inline void CSX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d half_v = _mm512_set1_pd(HALF);
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el2_real,el2_imag),_mm512_add_pd(el3_real,el3_imag)));
    __m512d sv_imag_pos2 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el2_imag,el2_real),_mm512_sub_pd(el3_imag,el3_real)));
    __m512d sv_real_pos3 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el2_real,el2_imag),_mm512_sub_pd(el3_real,el3_imag)));
    __m512d sv_imag_pos3 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el2_imag,el2_real),_mm512_add_pd(el3_imag,el3_real)));
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCSX Gate ================
inline void ConjugateCSX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d half_v = _mm512_set1_pd(HALF);
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el2_real,el2_imag),_mm512_sub_pd(el3_real,el3_imag)));
    __m512d sv_imag_pos2 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el2_imag,el2_real),_mm512_add_pd(el3_imag,el3_real)));
    __m512d sv_real_pos3 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_sub_pd(el2_real,el2_imag),_mm512_add_pd(el3_real,el3_imag)));
    __m512d sv_imag_pos3 = _mm512_mul_pd(half_v, _mm512_add_pd(_mm512_add_pd(el2_imag,el2_real),_mm512_sub_pd(el3_imag,el3_real)));
    STORE_Q1;
    OP_TAIL;
}
//============== CP Gate ================
inline void CP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta,
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d e3_real = _mm512_set1_pd(cos(theta));
    const __m512d e3_imag = _mm512_set1_pd(sin(theta));
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos1);
    const __m512d el3_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos3 = _mm512_sub_pd(_mm512_mul_pd(e3_real, el3_real), _mm512_mul_pd(e3_imag, el3_imag));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(e3_real, el3_imag), _mm512_mul_pd(e3_imag, el3_real));
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== ConjugateCP Gate ================
inline void ConjugateCP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    const __m512d e3_real = _mm512_set1_pd(cos(theta));
    const __m512d e3_imag = _mm512_set1_pd(-sin(theta));
    OP_HEAD_2Q;
    const __m512d el3_real = GET(sv_real, pos1);
    const __m512d el3_imag = GET(sv_imag, pos1);
    __m512d sv_real_pos3 = _mm512_sub_pd(_mm512_mul_pd(e3_real, el3_real), _mm512_mul_pd(e3_imag, el3_imag));
    __m512d sv_imag_pos3 = _mm512_add_pd(_mm512_mul_pd(e3_real, el3_imag), _mm512_mul_pd(e3_imag, el3_real));
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}

//============== CU Gate ================
inline void CU_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const ValType phi, const ValType lambda, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e1_real = _mm512_set1_pd(-cos(lambda)*sin(HALF*theta));
    const __m512d e1_imag = _mm512_set1_pd(-sin(lambda)*sin(HALF*theta));
    const __m512d e2_real = _mm512_set1_pd(cos(phi)*sin(HALF*theta));
    const __m512d e2_imag = _mm512_set1_pd(sin(phi)*sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(phi+lambda)*cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(sin(phi+lambda)*cos(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el2_real), _mm512_mul_pd(e1_real, el3_real)), 
                                          _mm512_mul_pd(e1_imag, el3_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el2_imag), _mm512_mul_pd(e1_real, el3_imag)), 
                                          _mm512_mul_pd(e1_imag, el3_real));
    __m512d sv_real_pos3 = _mm512_add_pd( _mm512_sub_pd(_mm512_mul_pd(e2_real, el2_real), _mm512_mul_pd(e2_imag, el2_imag)), 
                                          _mm512_sub_pd(_mm512_mul_pd(e3_real, el3_real), _mm512_mul_pd(e3_imag, el3_imag))); 
    __m512d sv_imag_pos3 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e2_real, el2_imag), _mm512_mul_pd(e2_imag, el2_real)), 
                                          _mm512_add_pd(_mm512_mul_pd(e3_real, el3_imag), _mm512_mul_pd(e3_imag, el3_real))); 
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCU Gate ================
inline void ConjugateCU_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const ValType phi, const ValType lambda, 
        const IdxType ctrl, const IdxType qubit)
{
    const __m512d e0_real = _mm512_set1_pd(cos(HALF*theta));
    const __m512d e1_real = _mm512_set1_pd(-cos(lambda)*sin(HALF*theta));
    const __m512d e1_imag = _mm512_set1_pd(sin(lambda)*sin(HALF*theta));
    const __m512d e2_real = _mm512_set1_pd(cos(phi)*sin(HALF*theta));
    const __m512d e2_imag = _mm512_set1_pd(-sin(phi)*sin(HALF*theta));
    const __m512d e3_real = _mm512_set1_pd(cos(phi+lambda)*cos(HALF*theta));
    const __m512d e3_imag = _mm512_set1_pd(-sin(phi+lambda)*cos(HALF*theta));
    OP_HEAD_2Q;
    LOAD_Q1;
    __m512d sv_real_pos2 = _mm512_sub_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el2_real), _mm512_mul_pd(e1_real, el3_real)), 
                                          _mm512_mul_pd(e1_imag, el3_imag));
    __m512d sv_imag_pos2 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e0_real, el2_imag), _mm512_mul_pd(e1_real, el3_imag)), 
                                          _mm512_mul_pd(e1_imag, el3_real));
    __m512d sv_real_pos3 = _mm512_add_pd( _mm512_sub_pd(_mm512_mul_pd(e2_real, el2_real), _mm512_mul_pd(e2_imag, el2_imag)), 
                                          _mm512_sub_pd(_mm512_mul_pd(e3_real, el3_real), _mm512_mul_pd(e3_imag, el3_imag))); 
    __m512d sv_imag_pos3 = _mm512_add_pd( _mm512_add_pd(_mm512_mul_pd(e2_real, el2_imag), _mm512_mul_pd(e2_imag, el2_real)), 
                                          _mm512_add_pd(_mm512_mul_pd(e3_real, el3_imag), _mm512_mul_pd(e3_imag, el3_real))); 
    STORE_Q1;
    OP_TAIL;
}
//============== ID Gate ================
inline void ID_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
}
//============== SWAP Gate ================
inline void SWAP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const __m512d el1_real = GET(sv_real, pos1);
    const __m512d el1_imag = GET(sv_imag, pos1);
    const __m512d el2_real = GET(sv_real, pos2);
    const __m512d el2_imag = GET(sv_imag, pos2);
    __m512d sv_real_pos1 = el2_real; 
    __m512d sv_imag_pos1 = el2_imag; 
    __m512d sv_real_pos2 = el1_real;
    __m512d sv_imag_pos2 = el1_imag;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    PUT(sv_real, pos2, sv_real_pos2);
    PUT(sv_imag, pos2, sv_imag_pos2);
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
    _Pragma("omp barrier")  
}

#endif
