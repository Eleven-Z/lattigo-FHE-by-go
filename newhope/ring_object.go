package newhope

import (
	"encoding/binary"
	"errors"
	"math/bits"
)

// Poly is the structure containing the coefficients of a polynomial.
type Poly struct {
	Coeffs [][]uint32 //Coefficients in CRT representation
}

// GetDegree returns the number of coefficients (degree) of the polynomial.
func (Pol *Poly) GetDegree() int {
	return len(Pol.Coeffs[0])
}

// GetLenModuli returns the number of modulies
func (Pol *Poly) GetLenModuli() int {
	return len(Pol.Coeffs)
}

// Zero sets all coefficient of the target polynomial to 0.
func (Pol *Poly) Zero() {
	for i := range Pol.Coeffs {
		for j := range Pol.Coeffs[0] {
			Pol.Coeffs[i][j] = 0
		}
	}
}

// CopyNew creates a new polynomial p1 which is a copy of the target polynomial.
func (Pol *Poly) CopyNew() (p1 *Poly) {
	p1 = new(Poly)
	p1.Coeffs = make([][]uint32, len(Pol.Coeffs))
	for i := range Pol.Coeffs {
		p1.Coeffs[i] = make([]uint32, len(Pol.Coeffs[i]))
		for j := range Pol.Coeffs[i] {
			p1.Coeffs[i][j] = Pol.Coeffs[i][j]
		}
	}

	return p1
}

// Copy copies the coefficients of p0 on p1 within the given context. Requiers p1 to be as big as the target context.
func (context *Context) Copy(p0, p1 *Poly) error {
	if len(p1.Coeffs) < len(context.Modulus) || uint32(len(p1.Coeffs[0])) < context.N {
		return errors.New("error : copy Poly, receiver poly is invalide")
	}
	for i := range context.Modulus {
		for j := uint32(0); j < context.N; j++ {
			p1.Coeffs[i][j] = p0.Coeffs[i][j]
		}
	}
	return nil
}

// Copy copies the coefficients of Pol on p1, require p1 to be at least as big as Pol.
func (Pol *Poly) Copy(p1 *Poly) error {
	if len(Pol.Coeffs) > len(p1.Coeffs) || len(Pol.Coeffs[0]) > len(p1.Coeffs[0]) {
		return errors.New("error : copy Poly, receiver poly is invalide")
	}
	for i := range Pol.Coeffs {
		for j := range Pol.Coeffs[i] {
			p1.Coeffs[i][j] = Pol.Coeffs[i][j]
		}
	}
	return nil
}

// WriteCoeffsTo converts a matrix of coefficients to a byte array.
func WriteCoeffsTo(pointer, N, numberModuli uint32, coeffs [][]uint32, data []byte) (uint32, error) {
	tmp := N << 3
	for i := uint32(0); i < numberModuli; i++ {
		for j := uint32(0); j < N; j++ {
			binary.BigEndian.PutUint32(data[pointer+(j<<3):pointer+((j+1)<<3)], coeffs[i][j])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffs converts a byte array to a matrix of coefficients.
func DecodeCoeffs(pointer, N, numberModuli uint32, coeffs [][]uint32, data []byte) (uint32, error) {
	tmp := N << 3
	for i := uint32(0); i < numberModuli; i++ {
		for j := uint32(0); j < N; j++ {
			coeffs[i][j] = binary.BigEndian.Uint32(data[pointer+(j<<3) : pointer+((j+1)<<3)])
		}
		pointer += tmp
	}

	return pointer, nil
}

// DecodeCoeffs converts a byte array to a matrix of coefficients.
func DecodeCoeffsNew(pointer, N, numberModuli uint32, coeffs [][]uint32, data []byte) (uint32, error) {
	tmp := N << 3
	for i := uint32(0); i < numberModuli; i++ {
		coeffs[i] = make([]uint32, N)
		for j := uint32(0); j < N; j++ {
			coeffs[i][j] = binary.BigEndian.Uint32(data[pointer+(j<<3) : pointer+((j+1)<<3)])
		}
		pointer += tmp
	}

	return pointer, nil
}

func (Pol *Poly) MarshalBinary() ([]byte, error) {

	N := uint32(len(Pol.Coeffs[0]))
	numberModulies := uint32(len(Pol.Coeffs))

	data := make([]byte, 2+((N*numberModulies)<<3))

	if numberModulies > 0xFF {
		return nil, errors.New("error : poly max modulies uint16 overflow")
	}

	data[0] = uint8(bits.Len32(uint32(N)) - 1)
	data[1] = uint8(numberModulies)

	var pointer uint32

	pointer = 2

	if _, err := WriteCoeffsTo(pointer, N, numberModulies, Pol.Coeffs, data); err != nil {
		return nil, err
	}

	return data, nil
}

func (Pol *Poly) UnMarshalBinary(data []byte) (*Poly, error) {

	N := uint32(int(1 << data[0]))
	numberModulies := uint32(int(data[1]))

	var pointer uint32

	pointer = 2

	if ((uint32(len(data)) - pointer) >> 3) != N*numberModulies {
		return nil, errors.New("error : invalid polynomial encoding")
	}

	if _, err := DecodeCoeffs(pointer, N, numberModulies, Pol.Coeffs, data); err != nil {
		return nil, err
	}

	return Pol, nil
}
