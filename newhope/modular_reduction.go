package newhope

import (
	"math/bits"
)

// https://eprint.iacr.org/2016/504.pdf
// Returns r = 3*C mod 12289
func kred(C int64) (r int64) {
	return 3* (C & 4095) - (C >> 12)
}

//============================
//=== MONTGOMERY REDUCTION ===
//============================

// mform returns a*2^64 mod q. It take thes input a in
// conventional form and return r which is the
// the montgomery form of a of a mod q with a radix of 2^32.
// Note : a * u can't overflow uint64 since a < q and u = floor(2^64-1/q)
func mform(a, q uint32, u uint64) (r uint32) {
	r = -uint32((uint64(a) * u)>>32) * q
	if r >= q {r -= q}
	return
}

// mformConstant is identical to mform, except that it runs in constant time
// and returns a value in [0, 2q-1] (it omits the conditional reduction).
// Note : a * u can't overflow uint64 since a < q and u = floor(2^64-1/q)
func mformconstant(a, q uint32, u uint64) (r uint32) {
	r = -uint32((uint64(a) * u)>>32) * q
	return
}

// invmform returns a*(1/2^32) mod q. It takes the input a in
// montgomery form mod q with a radix of 2^32 and returns r which is the normal form of a mod q.
func invmform(a, q, qinv uint32) (r uint32) {
	r, _ = bits.Mul32(a * qinv, q)
	r = q - r
	if r >= q {
		r -= q
	}
	return
}

// invmformConstant is indentical to invmform, except that it runs in constant time
// and returns a value in [0, 2q-1].
func invmformconstant(a, q, qinv uint32) (r uint32) {
	r, _ = bits.Mul32(a * qinv, q)
	r = q - r
	return
}

// mredparams computes the parameter qinv = (q^-1) mod 2^32,
// required for mred.
func mredparam(q uint32) (qinv uint32) {
	var x uint32
	qinv = 1
	x = q
	for i := 0; i < 31; i++ {
		qinv *= x
		qinv &= 0xFFFFFFFF
		x *= x
		x &= 0xFFFFFFFF
	}
	return
}

// mred computes x * y * (1/2^32) mod q. Requires that at least one of the inputs is in
// montgomery form. If only one of the inputs is in montgomery form (ex : a pre-computed constant),
// the result will be in normal form. If both inputs are in montgomery form, then the result
// will be in montgomery form.
func mred(x, y, q, qinv uint32) (r uint32) {
	ahi, alo := bits.Mul32(x, y)
	R := alo * qinv
	H, _ := bits.Mul32(R, q)
	r = ahi - H + q
	if r >= q {
		r -= q
	}
	return
}

// mredConstant is identical to mred except it runs in
// constant time and returns a value in [0, 2q-1].
func mredconstant(x, y, q, qinv uint32) (r uint32) {
	ahi, alo := bits.Mul32(x, y)
	R := alo * qinv
	H, _ := bits.Mul32(R, q)
	r = ahi - H + q
	return
}

//==========================
//=== BARRETT REDUCTION  ===
//==========================

// bredparams computes the parameters required for the bred with
// a radix of 2^64.
func bredparam(q uint32) (params uint64) {
	return 0xFFFFFFFFFFFFFFFF / uint64(q)
}

// bredadd reduces a 32 bit integer by q.
// Assumes that x <= 32bits. Useful when several additions
// are performed before a modular reduction, as it is faster than
// applying a conditional reduction after each addition.
func bredadd(x, q uint32, u uint64) (r uint32) {
	s0, _ := bits.Mul32(x, uint32(u>>32))	
	r = x - s0 * q
	if r >= q {r -= q}
	return
}

// bredaddconstant is indentical to BReAdd, except it runs
// in constant time and returns a value in [0, 2q-1].
func bredaddconstant(x, q uint32, u []uint32) uint32 {
	s0, _ := bits.Mul32(x, u[0])
	return x - s0*q
}

// bred compute a*b mod q for arbitrary a,b uint32. To be used
// when both a,b can not be pre-computed. However applying a montgomery
// transform on either a or b might be faster depending on the computation
// to do, especially if either a or b need to be multiplied with several other
// values.
func bred(x, y, q uint32, u uint64) (r uint32) {
	a := uint64(x) * uint64(y)
	m, _ := bits.Mul64(a, u)
	r = uint32(a - uint64(q) * m)
	if r >= q {r -= q}
	return
}

// bred is indentical to bred, except it runs
// in constant time and returns a value in [0, 2q-1].
func bredconstant(x, y, q uint32, u uint64) (r uint32) {
	a := uint64(x) * uint64(y)
	m, _ := bits.Mul64(a, u)
	r = uint32(a - uint64(q) * m)
	return
}

//===============================
//==== CONDITIONAL REDUCTION ====
//===============================

// cred reduce returns a mod q, where
// a is required to be in the range [0, 2q-1].
func cred(a, q uint32) uint32 {
	if a >= q {
		return a - q
	}
	return a
}
