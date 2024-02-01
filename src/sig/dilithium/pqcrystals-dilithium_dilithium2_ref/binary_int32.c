#include "binary_int32.h"
#include "rounding.h"
#include "params.h"

void int_masking_boolean_refresh(int32_t*pC, int32_t *pA,pMsm pMsx)
{
	for(int32_t j = 0; j < shares; ++j)
	{
		pC[j] = pA[j];
	}

	for(int32_t i = 0; i < shares;++i)
	{
		for(int32_t j = i+1; j < shares;++j)
		{
			int32_t r;
            pMsx.sm[i][j] = random_int32_t();
            r = pMsx.sm[i][j];
			pC[i] ^=r;
			pC[j] ^=r;
		}
	}
}

void int_masking_boolean_fullrefresh(int32_t*pC, int32_t *pA)
{
	for(int32_t j = 0; j < shares; ++j)
	{
		pC[j] = pA[j];
	}

	for(int32_t i = 0; i < shares;++i)
	{
		for(int32_t j = 1; j < shares;++j)
		{
			int32_t r;
			r = random_int32_t();
			pC[0] ^=r;
			pC[j] ^=r;
		}
	}
}

int32_t int_masking_boolean_fullxor(int32_t *pA)
{
	int32_t r ;
	int32_t x [shares];
	int_masking_boolean_fullrefresh(x,pA);
	r = x[0];
	for(int32_t j = 1; j < shares; ++j)
	{
		r ^=x[j];
	}
	return r;
}




void int_masking_boolean_xor(int32_t*pC, int32_t *pA, int32_t *pB)
{
	for(int32_t i = 0; i < shares; ++i)
	{
		pC[i] = (pA[i]) ^ (pB[i]);
	}
}

void int_masking_boolean_neg(int32_t*pC, int32_t *pA)
{
	for(int32_t i = 0; i < shares; ++i)
	{
		pC[i] = -pA[i];
	}
}

void int_masking_boolean_not(int32_t*pC, int32_t *pA)
{
	pC[0] = ~pA[0];
	for(int32_t i = 1; i < shares; ++i)
	{
		pC[i] = pA[i];
	}
}

void int_masking_boolean_lshift(int32_t*pC, int32_t *pA, int32_t shift)
{
	for(int32_t i = 0; i < shares; ++i)
	{
		pC[i] = pA[i] << shift;
	}
}

void int_masking_boolean_rshift(int32_t*pC, int32_t *pA, int32_t shift)
{
	for(int32_t i = 0; i < shares; ++i)
	{
		pC[i] = pA[i] >> shift;
	}
}






void int_masking_boolean_unprotected_generate(int32_t* pC,int32_t a)
{
	pC[0] = a;
	for(int32_t i = 1; i < shares;++i)
	{
		pC[i] = 0;
	}
}

int32_t int_masking_boolean_unprotected_recombine(int32_t* pA)
{
	int32_t r = pA[0];
	for (int32_t i = 1; i < shares; i++) r ^= pA[i];
	return r;
}

void int_masking_arithmetic_add(int32_t* pC,int32_t* pA,int32_t* pB)
{
	for(int32_t i = 0; i < shares; ++i) pC[i] = pA[i] + pB[i];
}

void int_masking_arithmetic_generate(int32_t* pC,int32_t a)
{
	pC[0] = a;
	for(int32_t i = 1; i < shares;++i)
	{
		pC[i] = random_int32_t();
		pC[0] -= pC[i];
	}
}

int32_t int_masking_arithmetic_unprotected_recombine(int32_t* pR)
{
	int32_t r = pR[0];
	for(int32_t i = 1; i < shares; ++i)
		r += pR[i];
	return r;
}


void int_masking_arithmetic_to_boolean_convert(int32_t*pC, int32_t *pA,pMtL*  pMR)
{
	int32_t b[shares];
	int_masking_boolean_generate(pC,pA[0],&pMR->pMm1);
	for(int32_t j = 1; j < shares;++j)
	{
		int_masking_boolean_generate( b,pA[j],&pMR->pMm2[j]);
		int_masking_boolean_add     (pC,pC,b ,&pMR->pMm3[j]);
	}
}


void int_masking_boolean_generate(int32_t* pC,int32_t a,svec *pVa)
{
    pC[0] = a;
    for(int32_t i = 1; i < shares;++i)
    {
        pVa->gen[i] = random_int32_t();
        pC[i] = pVa->gen[i];
        pC[0] ^= pC[i];
    }
}


void int_masking_boolean_add(int32_t*pC, int32_t *pA, int32_t *pB,pMtM *pMx)
{
    int32_t aux [shares];
    int32_t aux0[shares];
    int32_t p   [shares];
    int32_t g   [shares];
    const int32_t W = 5; /// log2 (w) = log2 (32)
    int_masking_boolean_xor(p,pA,pB);
    int_masking_boolean_and(g,pA,pB,pMx->pMs1);

    for(int32_t j = 1; j < W; ++j)
    {
        int32_t pow = 1 << (j-1);
        int_masking_boolean_lshift ( aux,   g, pow);
        int_masking_boolean_and    ( aux,   p, aux,pMx->pMs2[j]);
        int_masking_boolean_xor    (   g,   g, aux);
        int_masking_boolean_lshift (aux0,   p, pow);
        int_masking_boolean_refresh(aux0,aux0     ,pMx->pMs3[j]);
        int_masking_boolean_and    (   p,   p,aux0,pMx->pMs4[j]);
    }
    int_masking_boolean_lshift(aux, g,1<<(W-1));
    int_masking_boolean_and   (aux, p,     aux,pMx->pMs5);
    int_masking_boolean_xor   (  g, g,     aux);
    int_masking_boolean_lshift(aux, g,       1);
    int_masking_boolean_xor   ( pC,pA,      pB);
    int_masking_boolean_xor   ( pC,pC,     aux);
}


void int_masking_boolean_and(int32_t*pC, int32_t *pA, int32_t *pB,pMsm pMsx)
{
    int32_t r   [shares];
    for(int32_t i = 0; i < shares; ++i)
    {
        r[i] = (pA[i]) & (pB[i]);
    }

    for(int32_t i = 0; i < shares;++i)
    {
        for(int32_t j = i+1; j < shares;++j)
        {
            int32_t zij,zji,tmp;
            pMsx.sm[i][j] = random_int32_t();
            zij = pMsx.sm[i][j];
            zji  = ((pA[i])&(pB[j]));
            zji ^=  zij;
            tmp = ((pA[j])&(pB[i]));
            zji  ^= tmp;
            r[i] ^= zij;
            r[j] ^= zji;
        }
    }
    for(int32_t i = 0; i < shares; ++i)
    {
        pC[i] = r[i];
    }
}


//here is unmasking


//b_add7line
void int_unmasking_boolean_refresh(int32_t*pC, int32_t *pA,pMsm pMsx)
{
    for(int32_t j = 0; j < shares; ++j)
    {
        pC[j] = pA[j];
    }

    for(int32_t i = 0; i < shares;++i)
    {
        for(int32_t j = i+1; j < shares;++j)
        {
            int32_t r;
            r = pMsx.sm[i][j];
            pC[i] ^=r;
            pC[j] ^=r;
        }
    }
}


//con1line,2line
void int_unmasking_boolean_generate(int32_t* pC,int32_t a,svec *pVa)
{
    pC[0] = a;
    for(int32_t i = 1; i < shares;++i)
    {
        pC[i] = pVa->gen[i];
        pC[0] ^= pC[i];
    }
}

//b_add1,5,11,13,14line
void int_unmasking_boolean_xor(int32_t*pC, int32_t *pA, int32_t *pB)
{
    for(int32_t i = 0; i < shares; ++i)
    {
        pC[i] = (pA[i]) ^ (pB[i]);
    }
}

//b_add3,6,9,12line
void int_unmasking_boolean_lshift(int32_t*pC, int32_t *pA, int32_t shift)
{
    for(int32_t i = 0; i < shares; ++i)
    {
        pC[i] = pA[i] << shift;
    }
}

//HB3line
void int_unmasking_boolean_rshift(int32_t*pC, int32_t *pA, int32_t shift)
{
    for(int32_t i = 0; i < shares; ++i)
    {
        pC[i] = pA[i] >> shift;
    }
}




//b_add2,4,8,10line
void int_unmasking_boolean_and(int32_t*pC, int32_t *pA, int32_t *pB,pMsm pMsx)
{
    int32_t r   [shares];
    for(int32_t i = 0; i < shares; ++i)
    {
        r[i] = (pA[i]) & (pB[i]);
    }

    for(int32_t i = 0; i < shares;++i)
    {
        for(int32_t j = i+1; j < shares;++j)
        {
            int32_t zij,zji,tmp;
            zij = pMsx.sm[i][j];
            zji  = ((pA[i])&(pB[j]));
            zji ^=  zij;
            tmp = ((pA[j])&(pB[i]));
            zji  ^= tmp;
            r[i] ^= zij;
            r[j] ^= zji;
        }
    }
    for(int32_t i = 0; i < shares; ++i)
    {
        pC[i] = r[i];
    }
}

//HB1line
void int_unmasking_arithmetic_add(int32_t* pC,int32_t* pA,int32_t* pB)
{
    for(int32_t i = 0; i < shares; ++i) pC[i] = pA[i] + pB[i];
}


//con3line
void int_unmasking_boolean_add(int32_t*pC, int32_t *pA, int32_t *pB,pMtM *pMx)
{
    int32_t aux [shares];
    int32_t aux0[shares];
    int32_t p   [shares];
    int32_t g   [shares];
    const int32_t W = 5; /// log2 (w) = log2 (32)
    int_unmasking_boolean_xor(p,pA,pB);
    int_unmasking_boolean_and(g,pA,pB,pMx->pMs1);
    for(int32_t j = 1; j < W; ++j)
    {
        int32_t pow = 1 << (j-1);
        int_unmasking_boolean_lshift ( aux,   g, pow);
        int_unmasking_boolean_and    ( aux,   p, aux,pMx->pMs2[j]);
        int_unmasking_boolean_xor    (   g,   g, aux);
        int_unmasking_boolean_lshift (aux0,   p, pow);
        int_unmasking_boolean_refresh(aux0,aux0     ,pMx->pMs3[j]);
        int_unmasking_boolean_and    (   p,   p,aux0,pMx->pMs4[j]);
    }
    int_unmasking_boolean_lshift(aux, g,1<<(W-1));
    int_unmasking_boolean_and   (aux, p,     aux,pMx->pMs5);
    int_unmasking_boolean_xor   (  g, g,     aux);
    int_unmasking_boolean_lshift(aux, g,       1);
    int_unmasking_boolean_xor   ( pC,pA,      pB);
    int_unmasking_boolean_xor   ( pC,pC,     aux);
}

//HB2line
void int_unmasking_arithmetic_to_boolean_convert(int32_t*pC, int32_t *pA,pMtL*  pMR)
{
    int32_t b[shares];
    int_unmasking_boolean_generate(pC,pA[0],&pMR->pMm1);
    for(int32_t j = 1; j < shares;++j)
    {
        int_unmasking_boolean_generate( b,pA[j],&pMR->pMm2[j]);
        int_unmasking_boolean_add     (pC,pC,b ,&pMR->pMm3[j]);
    }
}










