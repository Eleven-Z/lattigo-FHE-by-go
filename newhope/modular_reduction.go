package newhope

import (
	"math/bits"
)

//============================
//=== MONTGOMERY REDUCTION ===
//============================

// MForm returns a*2^64 mod q. It take thes input a in
// conventional form and return r which is the
// the montgomery form of a of a mod q with a radix of 2^64.
func MForm(a, q uint32, u []uint32) (r uint32) {
	mhi, _ := bits.Mul32(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

// MFormConstant is identical to MForm, except that it runs in constant time
// and returns a value in [0, 2q-1] (it omits the conditional reduction).
func MFormConstant(a, q uint32, u []uint32) (r uint32) {
	mhi, _ := bits.Mul32(a, u[1])
	r = -(a*u[0] + mhi) * q
	return
}

// InvMForm returns a*(1/2^64) mod q. It takes the input a in
// montgomery form mod q with a radix of 2^ 64and returns r which is the normal form of a mod q.
func InvMForm(a, q, qInv uint32) (r uint32) {
	r, _ = bits.Mul32(a*qInv, q)
	r = q - r
	if r >= q {
		r -= q
	}
	return
}

// InvMFormConstant is indentical to InvMForm, except that it runs in constant time
// and returns a value in [0, 2q-1].
func InvMFormConstant(a, q, qInv uint32) (r uint32) {
	r, _ = bits.Mul32(a*qInv, q)
	r = q - r
	return
}

// MRedParams computes the parameter qInv = (q^-1) mod 2^32,
// required for MRed.
func MRedParams(q uint32) (qInv uint32) {
	var x uint32
	qInv = 1
	x = q
	for i := 0; i < 31; i++ {
		qInv *= x
		qInv &= 0xFFFFFFFF
		x *= x
		x &= 0xFFFFFFFF
	}
	return
}

// MRed computes x * y * (1/2^64) mod q. Requires that at least one of the inputs is in
// montgomery form. If only one of the inputs is in montgomery form (ex : a pre-computed constant),
// the result will be in normal form. If both inputs are in montgomery form, then the result
// will be in montgomery form.
func MRed(x, y, q, qInv uint32) (r uint32) {
	ahi, alo := bits.Mul32(x, y)
	R := alo * qInv
	H, _ := bits.Mul32(R, q)
	r = ahi - H + q
	if r >= q {
		r -= q
	}
	return
}

// MRedConstant is identical to MRed except it runs in
// constant time and returns a value in [0, 2q-1].
func MRedConstant(x, y, q, qInv uint32) (r uint32) {
	ahi, alo := bits.Mul32(x, y)
	R := alo * qInv
	H, _ := bits.Mul32(R, q)
	r = ahi - H + q
	return
}

//==========================
//=== BARRETT REDUCTION  ===
//==========================

// BRedParams computes the parameters required for the BRed with
// a radix of 2^64.
func BRedParams(q uint32) (params []uint32) {
	var bigR uint64
	bigR = 0xFFFFFFFFFFFFFFFF
	bigR /= uint64(q)

	// 2^radix // q
	mhi := uint32(bigR>>32)
	mlo := uint32(bigR)

	return []uint32{mhi, mlo}
}

// BRedAdd reduces a 64 bit integer by q.
// Assumes that x <= 64bits. Useful when several additions
// are performed before a modular reduction, as it is faster than
// applying a conditional reduction after each addition.
func BRedAdd(x, q uint32, u []uint32) (r uint32) {
	s0, _ := bits.Mul32(x, u[0])
	r = x - s0*q
	if r >= q {
		r -= q
	}
	return
}

// BRedAddConstant is indentical to BReAdd, except it runs
// in constant time and returns a value in [0, 2q-1].
func BRedAddConstant(x, q uint32, u []uint32) uint32 {
	s0, _ := bits.Mul32(x, u[0])
	return x - s0*q
}

// BRed compute a*b mod q for arbitrary a,b uint32. To be used
// when both a,b can not be pre-computed. However applying a montgomery
// transform on either a or b might be faster depending on the computation
// to do, especially if either a or b need to be multiplied with several other
// values.
func BRed(x, y, q uint32, u []uint32) (r uint32) {

	var lhi, mhi, mlo, s0, s1, carry uint32

	ahi, alo := bits.Mul32(x, y)

	// (alo*ulo)>>64

	lhi, _ = bits.Mul32(alo, u[1])

	// ((ahi*ulo + alo*uhi) + (alo*ulo))>>64

	mhi, mlo = bits.Mul32(alo, u[0])

	s0, carry = bits.Add32(mlo, lhi, 0)

	s1 = mhi + carry

	mhi, mlo = bits.Mul32(ahi, u[1])

	_, carry = bits.Add32(mlo, s0, 0)

	lhi = mhi + carry

	// (ahi*uhi) + (((ahi*ulo + alo*uhi) + (alo*ulo))>>64)

	s0 = ahi*u[0] + s1 + lhi

	r = alo - s0*q

	if r >= q {
		r -= q
	}

	return
}

// BRedConstant is indentical to BRed, except it runs
// in constant time and returns a value in [0, 2q-1].
func BRedConstant(x, y, q uint32, u []uint32) (r uint32) {

	var lhi, mhi, mlo, s0, s1, carry uint32

	ahi, alo := bits.Mul32(x, y)

	// alo*ulo

	lhi, _ = bits.Mul32(alo, u[1])

	// ahi*ulo + alo*uhi

	mhi, mlo = bits.Mul32(alo, u[0])

	s0, carry = bits.Add32(mlo, lhi, 0)

	s1 = mhi + carry

	mhi, mlo = bits.Mul32(ahi, u[1])

	_, carry = bits.Add32(mlo, s0, 0)

	lhi = mhi + carry

	// ahi*uhi

	s0 = ahi*u[0] + s1 + lhi

	r = alo - s0*q

	return
}

//===============================
//==== CONDITIONAL REDUCTION ====
//===============================

// CRed reduce returns a mod q, where
// a is required to be in the range [0, 2q-1].
func CRed(a, q uint32) uint32 {
	if a >= q {
		return a - q
	}
	return a
}
