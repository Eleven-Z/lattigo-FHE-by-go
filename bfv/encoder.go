package bfv

import (
	"errors"
	"github.com/ldsec/lattigo/ring"
	"math/bits"
)

// encoder is a structure storing the parameters encode values on a plaintext in a SIMD fashion.
type Encoder struct {
	indexMatrix  []uint64
	bfvcontext   *BfvContext
	simplescaler *ring.SimpleScaler
	polypool     *ring.Poly
}

// Newencoder creates a new encoder from the target bfvcontext.
func (bfvcontext *BfvContext) NewEncoder() (encoder *Encoder) {

	if bfvcontext.contextT.AllowsNTT() != true {
		panic("cannot create batch encoder : plaintext modulus does not allow NTT")
	}

	var m, gen, pos, index1, index2 uint64

	encoder = new(Encoder)

	encoder.bfvcontext = bfvcontext

	slots := bfvcontext.n

	encoder.indexMatrix = make([]uint64, slots)

	logN := uint64(bits.Len64(bfvcontext.n) - 1)

	row_size := bfvcontext.n >> 1
	m = (bfvcontext.n << 1)
	gen = bfvcontext.gen
	pos = 1

	for i := uint64(0); i < row_size; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		encoder.indexMatrix[i] = bitReverse64(index1, logN)
		encoder.indexMatrix[i|row_size] = bitReverse64(index2, logN)

		pos *= gen
		pos &= (m - 1)
	}

	encoder.simplescaler = ring.NewSimpleScaler(bfvcontext.t, bfvcontext.contextQ)
	encoder.polypool = bfvcontext.contextT.NewPoly()

	return encoder
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *Encoder) EncodeUint(coeffs []uint64, plaintext *Plaintext) error {

	if len(coeffs) > len(encoder.indexMatrix) {
		return errors.New("invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		return errors.New("invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder")
	}

	for i := 0; i < len(coeffs); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	if err := plaintext.EMB(encoder.bfvcontext); err != nil {
		return err
	}

	plaintext.Lift(encoder.bfvcontext)

	return nil
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. Also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient do not exceed half of the plaintext modulus.
func (encoder *Encoder) EncodeInt(coeffs []int64, plaintext *Plaintext) error {

	if len(coeffs) > len(encoder.indexMatrix) {
		return errors.New("invalid input to encode (number of coefficients must be smaller or equal to the context)")
	}

	if len(plaintext.value.Coeffs[0]) != len(encoder.indexMatrix) {
		return errors.New("invalid plaintext to receive encoding (number of coefficients does not match the context of the encoder)")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.bfvcontext.t) + coeffs[i])
		} else {
			plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		plaintext.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	if err := plaintext.EMB(encoder.bfvcontext); err != nil {
		return err
	}

	plaintext.Lift(encoder.bfvcontext)

	return nil
}

// DecodeUint decodes a batched plaintext and returns the coefficients in a uint64 slice.
func (encoder *Encoder) DecodeUint(plaintext *Plaintext) (coeffs []uint64) {

	encoder.simplescaler.Scale(plaintext.value, encoder.polypool)

	encoder.bfvcontext.contextT.NTT(encoder.polypool, encoder.polypool)

	coeffs = make([]uint64, encoder.bfvcontext.n)

	for i := uint64(0); i < encoder.bfvcontext.n; i++ {
		coeffs[i] = encoder.polypool.Coeffs[0][encoder.indexMatrix[i]]
	}

	return

}

// DecodeInt decodes a batched plaintext and returns the coefficients in an int64 slice. Also decodes the sign (by centering the values around the plaintext
// modulus).
func (encoder *Encoder) DecodeInt(plaintext *Plaintext) (coeffs []int64) {

	var value int64

	encoder.simplescaler.Scale(plaintext.value, encoder.polypool)

	encoder.bfvcontext.contextT.NTT(encoder.polypool, encoder.polypool)

	coeffs = make([]int64, encoder.bfvcontext.n)

	modulus := int64(encoder.bfvcontext.t)

	for i := uint64(0); i < encoder.bfvcontext.n; i++ {

		value = int64(encoder.polypool.Coeffs[0][encoder.indexMatrix[i]])

		coeffs[i] = value

		if value > modulus>>1 {
			coeffs[i] -= modulus
		}
	}

	return coeffs
}
