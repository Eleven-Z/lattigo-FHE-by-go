package newhope

import (
	"math/bits"
)

// Add adds p1 to p2 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Add(p1, p2, p3 *Poly) {
	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = cred(p1.Coeffs[i]+p2.Coeffs[i], context.q)
	}
}

// AddNoMod adds p1 to p2 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*context.q -1].
func (context *Context) AddNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = p1.Coeffs[i] + p2.Coeffs[i]
	}

}

// Sub subtract p2 to p1 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Sub(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = cred(p1.Coeffs[i]+(context.q-p2.Coeffs[i]), context.q)
	}
}

// SubNoMod subtract p2 to p1 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*context.q -1].
func (context *Context) SubNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = p1.Coeffs[i] + (context.q - p2.Coeffs[i])
	}
}

// Neg set all coefficient of p1 to there additive inverse, returning the result on p2.
func (context *Context) Neg(p1, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = context.q - p1.Coeffs[i]
	}

}

// Reduce applies a modular reduction over the coefficients of p1 returning the result on p2.
func (context *Context) Reduce(p1, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = bredadd(p1.Coeffs[i], context.q, context.bredparam)
	}
}

// Mod applies a modular reduction by m over the coefficients of p1, returning the result on p2.
func (context *Context) Mod(p1 *Poly, m uint32, p2 *Poly) {
	params := bredparam(m)

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = bredadd(p1.Coeffs[i], m, params)
	}
}

// AND applies a logical AND of m to the coefficients of p1,  returning the result on p2.
func (context *Context) AND(p1 *Poly, m uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = p1.Coeffs[i] & m
	}
}

// OR applies a logical OR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) OR(p1 *Poly, m uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = p1.Coeffs[i] | m
	}
}

// XOR applies a logical XOR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) XOR(p1 *Poly, m uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = p1.Coeffs[i] ^ m
	}
}

// MulCoeffs multiplies p1 by p2 coefficient wise with a Barrett modular reduction, returning the result on p3.
func (context *Context) MulCoeffs(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = bred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.bredparam)
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient wise with a Barret modular reduction, adding the result to p3 with modular reduction.
func (context *Context) MulCoeffsAndAdd(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = cred(p3.Coeffs[i]+bred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.bredparam), context.q)
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient wise with a Barrett modular reduction, adding the result to p3 without modular reduction.
func (context *Context) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] += bred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.bredparam)
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient wise with a montgomery modular reduction, returning the result on p3.
// Expects p1 and/or p2 to be in montgomery form for correctness (see mred).
func (context *Context) MulCoeffsMontgomery(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = mred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.mredparam)
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient wise with a montgomery modular reduction, adding the result to p3.
// Expects p1 and/or p2 to be in montgomery form for correctness (see mred).
func (context *Context) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = cred(p3.Coeffs[i]+mred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.mredparam), context.q)
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient wise with a montgomery modular reduction, adding the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see mred).
func (context *Context) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] += mred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.mredparam)
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient wise with a montgomery modular reduction, subtracting the result to p3 with modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see mred).
func (context *Context) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = cred(p3.Coeffs[i]+(context.q-mred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.mredparam)), context.q)
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient wise with a montgomery modular reduction, subtracting the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see mred).
func (context *Context) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = p3.Coeffs[i] + (context.q - mred(p1.Coeffs[i], p2.Coeffs[i], context.q, context.mredparam))
	}
}

// MulcoeffsConstant multiplies p1 by p2 coefficient wise with a constant time Barrett modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*context.q -1].
func (context *Context) MulCoeffsConstant(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = bredConstant(p1.Coeffs[i], p2.Coeffs[i], context.q, context.bredparam)
	}
}

// MulCoeffsConstantMontgomery multiplies p1 by p2 coefficient wise with a constant time Montgomery modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*context.q -1].
func (context *Context) MulCoeffsConstantMontgomery(p1, p2, p3 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p3.Coeffs[i] = mredconstant(p1.Coeffs[i], p2.Coeffs[i], context.q, context.mredparam)
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

	for i := uint32(0); i < context.n; i++ {

		for j := uint32(0); j < i; j++ {
			p3.Coeffs[j] = cred(p3.Coeffs[j]+(context.q-mred(p1Copy.Coeffs[i], p2Copy.Coeffs[context.n-i+j], context.q, context.mredparam)), context.q)
		}

		for j := uint32(i); j < context.n; j++ {
			p3.Coeffs[j] = cred(p3.Coeffs[j]+mred(p1Copy.Coeffs[i], p2Copy.Coeffs[j-i], context.q, context.mredparam), context.q)
		}
	}
}

// MulPolyNaiveMontgomery multiplies p1 by p2 with a naive convolution, returning the result on p3.
// Much faster than MulPolyNaive.
func (context *Context) MulPolyNaiveMontgomery(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	context.AND(p3, 0, p3)

	for i := uint32(0); i < context.n; i++ {

		for j := uint32(0); j < i; j++ {
			p3.Coeffs[j] = cred(p3.Coeffs[j]+(context.q-mred(p1Copy.Coeffs[i], p2Copy.Coeffs[context.n-i+j], context.q, context.mredparam)), context.q)
		}

		for j := uint32(i); j < context.n; j++ {
			p3.Coeffs[j] = cred(p3.Coeffs[j]+mred(p1Copy.Coeffs[i], p2Copy.Coeffs[j-i], context.q, context.mredparam), context.q)
		}
	}
}

// AddScalar adds to each coefficients of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) AddScalar(p1 *Poly, scalar uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = cred(p1.Coeffs[i]+scalar, context.q)
	}
}

// SubScalar subtracts to each coefficients of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) SubScalar(p1 *Poly, scalar uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = cred(p1.Coeffs[i]+(context.q-scalar), context.q)
	}
}

// MulScalar multiplies each coefficients of p1 by a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) MulScalar(p1 *Poly, scalar uint32, p2 *Poly) {
	var scalarMont uint32

	scalarMont = mform(bredadd(scalar, context.q, context.bredparam), context.q, context.bredparam)
	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = mred(p1.Coeffs[i], scalarMont, context.q, context.mredparam)
	}
}

// MForm sets p1 in conventional form to its montgomeryform, returning the result on p2.
func (context *Context) MForm(p1, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = mform(p1.Coeffs[i], context.q, context.bredparam)
	}
}

// MForm sets p1 in montgomeryform to its conventional form, returning the result on p2.
func (context *Context) InvMForm(p1, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p1.Coeffs[i] = invmform(p1.Coeffs[i], context.q, context.mredparam)
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

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = powerof2(p1.Coeffs[i], pow2, context.q, context.mredparam)
	}
}

// MulByVector multiplies p1 by a vector of uint32 coefficients and returns the result on p2.
func (context *Context) MulByVectorMontgomery(p1 *Poly, vector []uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] = mred(p1.Coeffs[i], vector[i], context.q, context.mredparam)
	}
}

// MulByVector multiplies p1 by a vector of uint32 coefficients and adds the result on p2 without modular reduction.
func (context *Context) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint32, p2 *Poly) {

	for i := uint32(0); i < context.n; i++ {
		p2.Coeffs[i] += mred(p1.Coeffs[i], vector[i], context.q, context.mredparam)
	}
}

// BitReverse applies a bit reverse permutation the coefficients of the input polynomial and returns the result on the receiver polynomial.
// Can safely be used for inplace permutation.
func (context *Context) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint32(bits.Len32(context.n) - 1)

	if p1 != p2 {

		for i := uint32(0); i < context.n; i++ {
			p2.Coeffs[bitreverse32(i, bitLenOfN)] = p1.Coeffs[i]
		}

	} else { // In place in case p1 = p2

		for i := uint32(0); i < context.n; i++ {
			j := bitreverse32(i, bitLenOfN)
			if i < j {
				p2.Coeffs[i], p2.Coeffs[j] = p2.Coeffs[j], p2.Coeffs[i]
			}
		}
	}
}
