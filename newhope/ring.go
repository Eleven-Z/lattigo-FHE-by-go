package newhope

import (
	"math/bits"
)

// Add adds p1 to p2 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Add(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p1.Coeffs[i][j]+p2.Coeffs[i][j], qi)
		}
	}
}

// AddNoMod adds p1 to p2 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) AddNoMod(p1, p2, p3 *Poly) {
	for i := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = p1.Coeffs[i][j] + p2.Coeffs[i][j]
		}
	}
}

// Sub subtract p2 to p1 coefficient wise and applies a modular reduction, returning the result on p3.
func (context *Context) Sub(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed((p1.Coeffs[i][j]+qi)-p2.Coeffs[i][j], qi)
		}
	}
}

// SubNoMod subtract p2 to p1 coefficient wise without modular reduction, returning the result on p3.
// The output range will be [0,2*Qi -1].
func (context *Context) SubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = (p1.Coeffs[i][j] + qi) - p2.Coeffs[i][j]
		}
	}
}

// Neg set all coefficient of p1 to there additive inverse, returning the result on p2.
func (context *Context) Neg(p1, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = qi - p1.Coeffs[i][j]
		}
	}
}

// Reduce applies a modular reduction over the coefficients of p1 returning the result on p2.
func (context *Context) Reduce(p1, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = BRedAdd(p1.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// Mod applies a modular reduction by m over the coefficients of p1, returning the result on p2.
func (context *Context) Mod(p1 *Poly, m uint32, p2 *Poly) {
	params := BRedParams(m)
	for i := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = BRedAdd(p1.Coeffs[i][j], m, params)
		}
	}
}

// AND applies a logical AND of m to the coefficients of p1,  returning the result on p2.
func (context *Context) AND(p1 *Poly, m uint32, p2 *Poly) {
	for i := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = p1.Coeffs[i][j] & m
		}
	}
}

// OR applies a logical OR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) OR(p1 *Poly, m uint32, p2 *Poly) {
	for i := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = p1.Coeffs[i][j] | m
		}
	}
}

// XOR applies a logical XOR of m to the coefficients of p1,  returning the result on p2.
func (context *Context) XOR(p1 *Poly, m uint32, p2 *Poly) {
	for i := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = p1.Coeffs[i][j] ^ m
		}
	}
}

// MulCoeffs multiplies p1 by p2 coefficient wise with a Barrett modular reduction, returning the result on p3.
func (context *Context) MulCoeffs(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = BRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient wise with a Barret modular reduction, adding the result to p3 with modular reduction.
func (context *Context) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p3.Coeffs[i][j]+BRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i]), qi)
		}
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient wise with a Barrett modular reduction, adding the result to p3 without modular reduction.
func (context *Context) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] += BRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient wise with a montgomery modular reduction, returning the result on p3.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])
		}
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient wise with a montgomery modular reduction, adding the result to p3.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p3.Coeffs[i][j]+MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i]), qi)
		}
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient wise with a montgomery modular reduction, adding the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] += MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])
		}
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient wise with a montgomery modular reduction, subtracting the result to p3 with modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = CRed(p3.Coeffs[i][j]+(qi-MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])), qi)
		}
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient wise with a montgomery modular reduction, subtracting the result to p3 without modular reduction.
// Expects p1 and/or p2 to be in montgomery form for correctness (see MRed).
func (context *Context) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = p3.Coeffs[i][j] + (qi - MRed(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i]))
		}
	}
}

// MulcoeffsConstant multiplies p1 by p2 coefficient wise with a constant time Barrett modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*Qi -1].
func (context *Context) MulCoeffsConstant(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = BRedConstant(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MulCoeffsConstantMontgomery multiplies p1 by p2 coefficient wise with a constant time Montgomery modular reduction, returning the result on p3.
// The output range of the modular reduction is [0, 2*Qi -1].
func (context *Context) MulCoeffsConstantMontgomery(p1, p2, p3 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p3.Coeffs[i][j] = MRedConstant(p1.Coeffs[i][j], p2.Coeffs[i][j], qi, context.mredParams[i])
		}
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

	for x, qi := range context.Modulus {

		for i := uint32(0); i < context.N; i++ {

			for j := uint32(0); j < i; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+(qi-MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][context.N-i+j], qi, context.mredParams[x])), qi)
			}

			for j := uint32(i); j < context.N; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][j-i], qi, context.mredParams[x]), qi)
			}
		}
	}
}

// MulPolyNaiveMontgomery multiplies p1 by p2 with a naive convolution, returning the result on p3.
// Much faster than MulPolyNaive.
func (context *Context) MulPolyNaiveMontgomery(p1, p2, p3 *Poly) {

	p1Copy := p1.CopyNew()
	p2Copy := p2.CopyNew()

	context.AND(p3, 0, p3)

	for x, qi := range context.Modulus {

		for i := uint32(0); i < context.N; i++ {

			for j := uint32(0); j < i; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+(qi-MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][context.N-i+j], qi, context.mredParams[x])), qi)
			}

			for j := uint32(i); j < context.N; j++ {
				p3.Coeffs[x][j] = CRed(p3.Coeffs[x][j]+MRed(p1Copy.Coeffs[x][i], p2Copy.Coeffs[x][j-i], qi, context.mredParams[x]), qi)
			}
		}
	}
}


// AddScalar adds to each coefficients of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) AddScalar(p1 *Poly, scalar uint32, p2 *Poly) {
	for i, Qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = CRed(p1.Coeffs[i][j]+scalar, Qi)
		}
	}
}

// SubScalar subtracts to each coefficients of p1 a scalar and applies a modular reduction, returing the result on p2.
func (context *Context) SubScalar(p1 *Poly, scalar uint32, p2 *Poly) {
	for i, Qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = CRed(p1.Coeffs[i][j]+(Qi-scalar), Qi)
		}
	}
}

// MulScalar multiplies each coefficients of p1 by a scalar and applies a modular reduction, returning the result on p2.
func (context *Context) MulScalar(p1 *Poly, scalar uint32, p2 *Poly) {
	var scalarMont uint32
	for i, Qi := range context.Modulus {
		scalarMont = MForm(BRedAdd(scalar, Qi, context.bredParams[i]), Qi, context.bredParams[i])
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = MRed(p1.Coeffs[i][j], scalarMont, Qi, context.mredParams[i])
		}
	}
}

// MForm sets p1 in conventional form to its montgomeryform, returning the result on p2.
func (context *Context) MForm(p1, p2 *Poly) {

	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = MForm(p1.Coeffs[i][j], qi, context.bredParams[i])
		}
	}
}

// MForm sets p1 in montgomeryform to its conventional form, returning the result on p2.
func (context *Context) InvMForm(p1, p2 *Poly) {

	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = InvMForm(p1.Coeffs[i][j], qi, context.mredParams[i])
		}
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
	for i, Qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = PowerOf2(p2.Coeffs[i][j], pow2, Qi, context.mredParams[i])
		}
	}
}



// MulByVector multiplies p1 by a vector of uint32 coefficients and returns the result on p2.
func (context *Context) MulByVectorMontgomery(p1 *Poly, vector []uint32, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] = MRed(p1.Coeffs[i][j], vector[j], qi, context.mredParams[i])
		}
	}
}

// MulByVector multiplies p1 by a vector of uint32 coefficients and adds the result on p2 without modular reduction.
func (context *Context) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint32, p2 *Poly) {
	for i, qi := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p2.Coeffs[i][j] += MRed(p1.Coeffs[i][j], vector[j], qi, context.mredParams[i])
		}
	}
}

// BitReverse applies a bit reverse permutation the coefficients of the input polynomial and returns the result on the receiver polynomial.
// Can safely be used for inplace permutation.
func (context *Context) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint32(bits.Len32(context.N) - 1)

	if p1 != p2 {
		for i := range context.Modulus {
			for j := uint32(0); j < context.N; j++ {
				p2.Coeffs[i][bitReverse32(j, bitLenOfN)] = p1.Coeffs[i][j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range context.Modulus {
			for i := uint32(0); i < context.N; i++ {
				j := bitReverse32(i, bitLenOfN)
				if i < j {
					p2.Coeffs[x][i], p2.Coeffs[x][j] = p2.Coeffs[x][j], p2.Coeffs[x][i]
				}
			}
		}
	}
}