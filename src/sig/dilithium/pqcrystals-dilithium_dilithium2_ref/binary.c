#include "binary.h"
#include "params.h"

void masking_boolean_refresh(uint32_t*pC, uint32_t *pA)
{
	for(uint32_t j = 0; j < shares; ++j)
	{
		pC[j] = pA[j];
	}

	for(uint32_t i = 0; i < shares;++i)
	{
		for(uint32_t j = i+1; j < shares;++j)
		{
			uint32_t r;
			r = rand_uint32_t();
			pC[i] ^=r;
			pC[j] ^=r;
		}
	}
}

void masking_boolean_fullrefresh(uint32_t*pC, uint32_t *pA)
{
	for(uint32_t j = 0; j < shares; ++j)
	{
		pC[j] = pA[j];
	}

	for(uint32_t i = 0; i < shares;++i)
	{
		for(uint32_t j = 1; j < shares;++j)
		{
			uint32_t r;
			r = rand_uint32_t();
			pC[0] ^=r;
			pC[j] ^=r;
		}
	}
}

uint32_t masking_boolean_fullxor(uint32_t *pA)
{
	uint32_t r ;
	uint32_t x [shares];
	masking_boolean_fullrefresh(x,pA);
	r = x[0];
	for(uint32_t j = 1; j < shares; ++j)
	{
		r ^=x[j];
	}
	return r;
}

void masking_boolean_generate(uint32_t* pC,uint32_t a)
{
	pC[0] = a;
	for(uint32_t i = 1; i < shares;++i)
	{
		pC[i] = rand_uint32_t();
		pC[0] ^= pC[i];
	}
}

void masking_boolean_xor(uint32_t*pC, uint32_t *pA, uint32_t *pB)
{
	for(uint32_t i = 0; i < shares; ++i)
	{
		pC[i] = (pA[i]) ^ (pB[i]);
	}
}

void masking_boolean_neg(uint32_t*pC, uint32_t *pA)
{
	for(uint32_t i = 0; i < shares; ++i)
	{
		pC[i] = -pA[i];
	}
}

void masking_boolean_not(uint32_t*pC, uint32_t *pA)
{
	pC[0] = ~pA[0];
	for(uint32_t i = 1; i < shares; ++i)
	{
		pC[i] = pA[i];
	}
}

void masking_boolean_lshift(uint32_t*pC, uint32_t *pA, uint32_t shift)
{
	for(uint32_t i = 0; i < shares; ++i)
	{
		pC[i] = pA[i] << shift;
	}
}

void masking_boolean_rshift(uint32_t*pC, uint32_t *pA, uint32_t shift)
{
	for(uint32_t i = 0; i < shares; ++i)
	{
		pC[i] = pA[i] >> shift;
	}
}

void masking_boolean_and(uint32_t*pC, uint32_t *pA, uint32_t *pB)
{
	uint32_t r   [shares];
	for(uint32_t i = 0; i < shares; ++i)
	{
		r[i] = (pA[i]) & (pB[i]);
	}

	for(uint32_t i = 0; i < shares;++i)
	{
		for(uint32_t j = i+1; j < shares;++j)
		{
			uint32_t zij,zji,tmp;
			zij = rand_uint32_t();
			zji  = ((pA[i])&(pB[j]));
			zji ^=  zij;
			tmp = ((pA[j])&(pB[i]));
			zji  ^= tmp;
			r[i] ^= zij;
			r[j] ^= zji;
		}
	}
	for(uint32_t i = 0; i < shares; ++i)
	{
		pC[i] = r[i];
	}
}

void masking_boolean_add(uint32_t*pC, uint32_t *pA, uint32_t *pB)
{
	uint32_t aux [shares];
	uint32_t aux0[shares];
	uint32_t p   [shares];
	uint32_t g   [shares];
	const uint32_t W = 5; /// log2 (w) = log2 (32)
	masking_boolean_xor(p,pA,pB);
	masking_boolean_and(g,pA,pB);
	for(uint32_t j = 1; j < W; ++j)
	{
		uint32_t pow = 1 << (j-1);
		masking_boolean_lshift ( aux,   g, pow);
		masking_boolean_and    ( aux,   p, aux);
		masking_boolean_xor    (   g,   g, aux);
		masking_boolean_lshift (aux0,   p, pow);
		masking_boolean_refresh(aux0,aux0     );
		masking_boolean_and    (   p,   p,aux0);
	}
	masking_boolean_lshift(aux, g,1<<(W-1));
	masking_boolean_and   (aux, p,     aux);
	masking_boolean_xor   (  g, g,     aux);
	masking_boolean_lshift(aux, g,       1);
	masking_boolean_xor   ( pC,pA,      pB);
	masking_boolean_xor   ( pC,pC,     aux);
}


void masking_boolean_unprotected_generate(uint32_t* pC,uint32_t a)
{
	pC[0] = a;
	for(uint32_t i = 1; i < shares;++i)
	{
		pC[i] = 0;
	}
}

uint32_t masking_boolean_unprotected_recombine(uint32_t* pA)
{
	uint32_t r = pA[0];
	for (uint32_t i = 1; i < shares; i++) r ^= pA[i];
	return r;
}

void masking_arithmetic_add(uint32_t* pC,uint32_t* pA,uint32_t* pB)
{
	for(uint32_t i = 0; i < shares; ++i) pC[i] = pA[i] + pB[i];
}

void masking_arithmetic_generate(uint32_t* pC,uint32_t a)
{
	pC[0] = a;
	for(uint32_t i = 1; i < shares;++i)
	{
		pC[i] = rand_uint32_t();
		pC[0] -= pC[i];
	}
}

uint32_t masking_arithmetic_unprotected_recombine(uint32_t* pR)
{
	uint32_t r = pR[0];
	for(uint32_t i = 1; i < shares; ++i)
		r += pR[i];
	return r;
}


void masking_arithmetic_to_boolean_convert(uint32_t*pC, uint32_t *pA)
{
	uint32_t b[shares];
	masking_boolean_generate(pC,pA[0]);
	for(uint32_t j = 1; j < shares;++j)
	{
		masking_boolean_generate( b,pA[j]);
		masking_boolean_add     (pC,pC,b );
	}
}

void masking_boolean_to_arithmetic_convert(uint32_t*pC, uint32_t *pA)
{
	uint32_t a[shares];
	uint32_t y[shares];
	uint32_t z[shares];
	for(uint32_t i = 0; i < shares-1;++i) pC[i] = rand_uint32_t();
	for(uint32_t i = 0; i < shares-1;++i) a[i] = -pC[i];
	a[shares-1] = 0;
	masking_arithmetic_to_boolean_convert(y,a);
	masking_boolean_add(z,pA,y);
	pC[shares-1] = masking_boolean_fullxor(z);
}
