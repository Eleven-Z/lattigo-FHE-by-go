package newhope

import (
	"math/bits"
)

// Add adds p1 to p2 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Add(p1, p2, p3 *Poly) {
	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = CRed(p1.Coeffs[i]+p2.Coeffs[i], context.Modulus)
	}
}

// AddNoMod adds p1 to p2 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*context.Modulus -1].
func (context *Context) AddNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = p1.Coeffs[i] + p2.Coeffs[i]
	}

}

// Sub subtract p2 to p1 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Sub(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = CRed(p1.Coeffs[i]+(context.Modulus-p2.Coeffs[i]), context.Modulus)
	}
}

// SubNoMod subtract p2 to p1 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*context.Modulus -1].
func (context *Context) SubNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = p1.Coeffs[i] + (context.Modulus - p2.Coeffs[i])
	}
}

// Neg set all coefficient of p1 to there additive inverse, returning the result on p2.
func (context *Context) Neg(p1, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = context.Modulus - p1.Coeffs[i]
	}

}

// Reduce applies a modular reduction over the coefficients of p1 returning the result on p2.
func (context *Context) Reduce(p1, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = BRedAdd(p1.Coeffs[i], context.Modulus, context.bredParams)
	}
}

// Mod applies a modular reduction by m over the coefficients of p1, returning the result on p2.
func (context *Context) Mod(p1 *Poly, m uint32, p2 *Poly) {
	params := BRedParams(m)

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = BRedAdd(p1.Coeffs[i], m, params)
	}
}

// AND applies a logical AND of m to the coefficients of p1,  returning the result on p2.
func (context *Context) AND(p1 *Poly, m uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = p1.Coeffs[i] & m
	}
}

// OR applies a logical OR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) OR(p1 *Poly, m uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = p1.Coeffs[i] | m
	}
}

// XOR applies a logical XOR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) XOR(p1 *Poly, m uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = p1.Coeffs[i] ^ m
	}
}

// MulCoeffs multiplies p1 by p2 coefficient wise with a Barrett modular reduction, returning the result on p3.
func (context *Context) MulCoeffs(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = BRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.bredParams)
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient wise with a Barret modular reduction, adding the result to p3 with modular reduction.
func (context *Context) MulCoeffsAndAdd(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = CRed(p3.Coeffs[i]+BRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.bredParams), context.Modulus)
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient wise with a Barrett modular reduction, adding the result to p3 without modular reduction.
func (context *Context) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] += BRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.bredParams)
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient wise with a montgomery modular reduction, returning the result on p3.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomery(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = MRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.mredParams)
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient wise with a montgomery modular reduction, adding the result to p3.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = CRed(p3.Coeffs[i]+MRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.mredParams), context.Modulus)
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient wise with a montgomery modular reduction, adding the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] += MRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.mredParams)
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient wise with a montgomery modular reduction, subtracting the result to p3 with modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = CRed(p3.Coeffs[i]+(context.Modulus-MRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.mredParams)), context.Modulus)
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient wise with a montgomery modular reduction, subtracting the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = p3.Coeffs[i] + (context.Modulus - MRed(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.mredParams))
	}
}

// MulcoeffsConstant multiplies p1 by p2 coefficient wise with a constant time Barrett modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*context.Modulus -1].
func (context *Context) MulCoeffsConstant(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = BRedConstant(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.bredParams)
	}
}

// MulCoeffsConstantMontgomery multiplies p1 by p2 coefficient wise with a constant time Montgomery modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*context.Modulus -1].
func (context *Context) MulCoeffsConstantMontgomery(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p3.Coeffs[i] = MRedConstant(p1.Coeffs[i], p2.Coeffs[i], context.Modulus, context.mredParams)
	}
}

// MulPoly multiplies p1 by p2 and returns the result on p3.
func (context *Context) MulPoly(p1, p2, p3 *Poly) {

	a := context.NewPoly()
	b := context.NewPoly()

	context.NTT(p1, a)
	context.NTT(p2, b)
	context.MulCoeffs(a, b, p3)
	context.InvNTT(p3, p3)
}

// MulPolyMontgomery multiplies p1 by p2 and returns the result on p3.
// Expect wither p1 or p2 to be in montgomery form for correctness.
func (context *Context) MulPolyMontgomery(p1, p2, p3 *Poly) {

	a := context.NewPoly()
	b := context.NewPoly()

	context.NTT(p1, a)
	context.NTT(p2, b)
	context.MulCoeffsMontgomery(a, b, p3)
	context.InvNTT(p3, p3)
}

// MulPolyNaive multiplies p1 by p2 with a naive convolution, returning the result on p3.
func (context *Context) MulPolyNaive(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	context.MForm(p1Copy, p1Copy)

	context.AND(p3, 0, p3)

	for i := uint32(0); i < context.N; i++ {

		for j := uint32(0); j < i; j++ {
			p3.Coeffs[j] = CRed(p3.Coeffs[j]+(context.Modulus-MRed(p1Copy.Coeffs[i], p2Copy.Coeffs[context.N-i+j], context.Modulus, context.mredParams)), context.Modulus)
		}

		for j := uint32(i); j < context.N; j++ {
			p3.Coeffs[j] = CRed(p3.Coeffs[j]+MRed(p1Copy.Coeffs[i], p2Copy.Coeffs[j-i], context.Modulus, context.mredParams), context.Modulus)
		}
	}
}

// MulPolyNaiveMontgomery multiplies p1 by p2 with a naive convolution, returning the result on p3.
// Much faster than MulPolyNaive.
func (context *Context) MulPolyNaiveMontgomery(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	context.AND(p3, 0, p3)

	for i := uint32(0); i < context.N; i++ {

		for j := uint32(0); j < i; j++ {
			p3.Coeffs[j] = CRed(p3.Coeffs[j]+(context.Modulus-MRed(p1Copy.Coeffs[i], p2Copy.Coeffs[context.N-i+j], context.Modulus, context.mredParams)), context.Modulus)
		}

		for j := uint32(i); j < context.N; j++ {
			p3.Coeffs[j] = CRed(p3.Coeffs[j]+MRed(p1Copy.Coeffs[i], p2Copy.Coeffs[j-i], context.Modulus, context.mredParams), context.Modulus)
		}
	}
}

// AddScalar adds to each coefficients of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) AddScalar(p1 *Poly, scalar uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = CRed(p1.Coeffs[i]+scalar, context.Modulus)
	}
}

// SubScalar subtracts to each coefficients of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) SubScalar(p1 *Poly, scalar uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = CRed(p1.Coeffs[i]+(context.Modulus-scalar), context.Modulus)
	}
}

// MulScalar multiplies each coefficients of p1 by a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) MulScalar(p1 *Poly, scalar uint32, p2 *Poly) {
	var scalarMont uint32

	scalarMont = MForm(BRedAdd(scalar, context.Modulus, context.bredParams), context.Modulus, context.bredParams)
	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = MRed(p1.Coeffs[i], scalarMont, context.Modulus, context.mredParams)
	}
}

// MForm sets p1 in conventional form to its montgomeryform, returning the result on p2.
func (context *Context) MForm(p1, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = MForm(p1.Coeffs[i], context.Modulus, context.bredParams)
	}
}

// MForm sets p1 in montgomeryform to its conventional form, returning the result on p2.
func (context *Context) InvMForm(p1, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p1.Coeffs[i] = InvMForm(p1.Coeffs[i], context.Modulus, context.mredParams)
	}
}

// MulByPow2New multiplies the input polynomial by 2^pow2 and returns the result on a new polynomial.
func (context *Context) MulByPow2New(p1 *Poly, pow2 uint32) (p2 *Poly) {
	p2 = context.NewPoly()
	context.MulByPow2(p1, pow2, p2)
	return
}

// MulByPow2 multiplies the input polynomial by 2^pow2 and returns the result on the receiver polynomial.
func (context *Context) MulByPow2(p1 *Poly, pow2 uint32, p2 *Poly) {
	context.MForm(p1, p2)

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = PowerOf2(p1.Coeffs[i], pow2, context.Modulus, context.mredParams)
	}
}

// MulByVector multiplies p1 by a vector of uint32 coefficients and returns the result on p2.
func (context *Context) MulByVectorMontgomery(p1 *Poly, vector []uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] = MRed(p1.Coeffs[i], vector[i], context.Modulus, context.mredParams)
	}
}

// MulByVector multiplies p1 by a vector of uint32 coefficients and adds the result on p2 without modular reduction.
func (context *Context) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint32, p2 *Poly) {

	for i := uint32(0); i < context.N; i++ {
		p2.Coeffs[i] += MRed(p1.Coeffs[i], vector[i], context.Modulus, context.mredParams)
	}
}

// BitReverse applies a bit reverse permutation the coefficients of the input polynomial and returns the result on the receiver polynomial.
// Can safely be used for inplace permutation.
func (context *Context) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint32(bits.Len32(context.N) - 1)

	if p1 != p2 {

		for i := uint32(0); i < context.N; i++ {
			p2.Coeffs[bitReverse32(i, bitLenOfN)] = p1.Coeffs[i]
		}

	} else { // In place in case p1 = p2

		for i := uint32(0); i < context.N; i++ {
			j := bitReverse32(i, bitLenOfN)
			if i < j {
				p2.Coeffs[i], p2.Coeffs[j] = p2.Coeffs[j], p2.Coeffs[i]
			}
		}
	}
}
