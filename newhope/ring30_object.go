package newhope

import (
	"encoding/binary"
	"errors"
	"math/bits"
)

// Poly is the structure containing the coefficients of a polynomial.
type Poly struct {
	Coeffs []uint32 //Coefficients in CRT representation
}

// GetDegree returns the number of coefficients (degree) of the polynomial.
func (Pol *Poly) GetDegree() int {
	return len(Pol.Coeffs)
}

// Zero sets all coefficient of the target polynomial to 0.
func (Pol *Poly) Zero() {
	for i := range Pol.Coeffs {
		Pol.Coeffs[i] = 0
	}
}

// CopyNew creates a new polynomial p1 which is a copy of the target polynomial.
func (Pol *Poly) CopyNew() (p1 *Poly) {
	p1 = new(Poly)
	p1.Coeffs = make([]uint32, len(Pol.Coeffs))
	for i := range Pol.Coeffs {
		p1.Coeffs[i] = Pol.Coeffs[i]

	}
	return p1
}

// Copy copies the coefficients of p0 on p1 within the given context. Requiers p1 to be as big as the target context.
func (context *Context) Copy(p0, p1 *Poly) error {
	if uint32(len(p1.Coeffs)) < context.N {
		return errors.New("error : copy Poly, receiver poly is invalide")
	}

	for i := uint32(0); i < context.N; i++ {
		p1.Coeffs[i] = p0.Coeffs[i]
	}

	return nil
}

// Copy copies the coefficients of Pol on p1, require p1 to be at least as big as Pol.
func (Pol *Poly) Copy(p1 *Poly) error {
	if len(Pol.Coeffs) > len(p1.Coeffs) {
		return errors.New("error : copy Poly, receiver poly is invalide")
	}

	for i := range Pol.Coeffs {
		p1.Coeffs[i] = Pol.Coeffs[i]
	}

	return nil
}

// WriteCoeffsTo converts a matrix of coefficients to a byte array.
func WriteCoeffsTo(pointer, N uint32, coeffs []uint32, data []byte) (uint32, error) {

	for i := uint32(0); i < N; i++ {
		binary.BigEndian.PutUint32(data[pointer+(i<<2):pointer+((i+1)<<2)], coeffs[i])
	}

	return pointer, nil
}

// DecodeCoeffs converts a byte array to a matrix of coefficients.
func DecodeCoeffs(pointer, N uint32, coeffs []uint32, data []byte) (uint32, error) {

	for i := uint32(0); i < N; i++ {
		coeffs[i] = binary.BigEndian.Uint32(data[pointer+(i<<2) : pointer+((i+1)<<2)])
	}

	return pointer, nil
}

func (Pol *Poly) MarshalBinary() ([]byte, error) {

	N := uint32(len(Pol.Coeffs))

	data := make([]byte, 1+(N<<2))

	data[0] = uint8(bits.Len32(uint32(N)) - 1)

	var pointer uint32

	pointer = 1

	if _, err := WriteCoeffsTo(pointer, N, Pol.Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

func (Pol *Poly) UnMarshalBinary(data []byte) (*Poly, error) {

	N := uint32(1 << data[0])

	var pointer uint32

	pointer = 1

	if ((uint32(len(data)) - pointer) >> 2) != N {
		return nil, errors.New("error : invalid polynomial encoding")
	}

	if _, err := DecodeCoeffs(pointer, N, Pol.Coeffs, data); err != nil {
		return nil, err
	}

	return Pol, nil
}
