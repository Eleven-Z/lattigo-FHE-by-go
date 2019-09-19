// Package ring implelents a RNS-accelerated modular arithmetic operations for polynomials, including: RNS basis extension; RNS rescaling;  number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package newhope

import (
	"bytes"
	"crypto/rand"
	"encoding/binary"
	"encoding/gob"
	"errors"
	"math/bits"
)

//==============================
//===== POLYNOMIAL CONTEXT =====
//==============================

// Context is a structure keeping all the variable required to operate on a polynomial represented in this context.
// This include its moduli, crt reconstruction, modular reduction and ntt transformation.
type Context struct {

	// Polynomial nb.Coefficients
	N uint32

	// Moduli
	Modulus []uint32

	// 2^bit_length(Qi) - 1
	mask []uint32

	// Determines if NTT can be used with the current context.
	validated bool

	// Fast reduction parameters
	bredParams [][]uint32
	mredParams []uint32

	//NTT Parameters
	psiMont    []uint32 //2nth primitive root in montgomery form
	psiInvMont []uint32 //2nth inverse primitive root in montgomery form

	nttPsi    [][]uint32 //powers of the inverse of the 2nth primitive root in montgomery form (in bitreversed order)
	nttPsiInv [][]uint32 //powers of the inverse of the 2nth primitive root in montgomery form (in bitreversed order)
	nttNInv   []uint32   //[N^-1] mod Qi in montgomery form
}

// NewContext generates a new empty context.
func NewContext() *Context {
	return new(Context)
}

// SetParameters initialize the parameters of an empty context with N and the provided moduli.
// Only checks that N is a power of 2 and computes all the variable that aren't used for the NTT.
func (context *Context) SetParameters(N uint32, Modulus []uint32) error {

	// Checks if N is a power of 2
	if (N&(N-1)) != 0 && N != 0 {
		return errors.New("invalid ring degree (must be a power of 2)")
	}

	context.validated = false

	context.N = N

	context.Modulus = make([]uint32, len(Modulus))
	context.mask = make([]uint32, len(Modulus))

	for i, qi := range Modulus {
		context.Modulus[i] = qi
		context.mask[i] = (1 << uint32(bits.Len32(qi))) - 1
	}

	// Computes the fast reduction parameters
	context.bredParams = make([][]uint32, len(context.Modulus))
	context.mredParams = make([]uint32, len(context.Modulus))

	for i, qi := range context.Modulus {

		//Computes the fast modular reduction parameters for the Context
		context.bredParams[i] = BRedParams(qi)

		// If qi is not a power of 2, we can compute the MRedParams (else it should not
		// because it will return an error and there is no valid montgomery form mod a power of 2)
		if (qi&(qi-1)) != 0 && qi != 0 {
			context.mredParams[i] = MRedParams(qi)
		}
	}

	return nil
}

// ValidateParameters checks that N has beed correctly initialized, and checks that each moduli is a prime congruent to 1 mod 2N (i.e. allowing NTT).
// Then it computes the variables required for the NTT. ValidateParameters purpose is to validate that the moduli allow the NTT and compute the
// NTT parameters.
func (context *Context) ValidateParameters() error {

	if context.validated {
		return nil
	}

	if context.N == 0 || context.Modulus == nil {
		return errors.New("error : invalid context parameters (missing)")
	}

	// CHECKS IF VALIDE NTT
	// Checks if each qi is Prime and if qi = 1 mod 2n
	for _, qi := range context.Modulus {
		if IsPrime(qi) == false || qi&((context.N<<1)-1) != 1 {
			context.validated = false
			return errors.New("warning : provided modulus does not allow NTT")
		}
	}

	context.psiMont = make([]uint32, len(context.Modulus))
	context.psiInvMont = make([]uint32, len(context.Modulus))
	context.nttPsi = make([][]uint32, len(context.Modulus))
	context.nttPsiInv = make([][]uint32, len(context.Modulus))
	context.nttNInv = make([]uint32, len(context.Modulus))

	bitLenofN := uint32(bits.Len32(context.N) - 1)

	for i, qi := range context.Modulus {

		//2.1 Computes N^(-1) mod Q in Montgomery form
		context.nttNInv[i] = MForm(ModExp(context.N, qi-2, qi), qi, context.bredParams[i])

		//2.2 Computes Psi and PsiInv in Montgomery form
		context.nttPsi[i] = make([]uint32, context.N)
		context.nttPsiInv[i] = make([]uint32, context.N)

		//Finds a 2nth primitive Root
		g := primitiveRoot(qi)

		_2n := uint32(context.N << 1)

		power := (qi - 1) / _2n
		powerInv := (qi - 1) - power

		//Computes Psi and PsiInv in Montgomery Form
		PsiMont := MForm(ModExp(g, power, qi), qi, context.bredParams[i])
		PsiInvMont := MForm(ModExp(g, powerInv, qi), qi, context.bredParams[i])

		context.psiMont[i] = PsiMont
		context.psiInvMont[i] = PsiInvMont

		context.nttPsi[i][0] = MForm(1, qi, context.bredParams[i])
		context.nttPsiInv[i][0] = MForm(1, qi, context.bredParams[i])

		// Computes nttPsi[j] = nttPsi[j-1]*Psi and nttPsiInv[j] = nttPsiInv[j-1]*PsiInv
		for j := uint32(1); j < context.N; j++ {

			indexReversePrev := bitReverse32(j-1, bitLenofN)
			indexReverseNext := bitReverse32(j, bitLenofN)

			context.nttPsi[i][indexReverseNext] = MRed(context.nttPsi[i][indexReversePrev], PsiMont, qi, context.mredParams[i])
			context.nttPsiInv[i][indexReverseNext] = MRed(context.nttPsiInv[i][indexReversePrev], PsiInvMont, qi, context.mredParams[i])
		}
	}

	context.validated = true

	return nil
}

// Used to export the context. Minimal information to recover the full context.
type smallContext struct {
	N       uint32
	Modulus []uint32
}

func (context *Context) MarshalBinary() ([]byte, error) {

	parameters := smallContext{context.N, context.Modulus}

	var buf bytes.Buffer
	enc := gob.NewEncoder(&buf)
	if err := enc.Encode(parameters); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

func (context *Context) UnMarshalBinary(data []byte) error {

	parameters := smallContext{}

	reader := bytes.NewReader(data)
	dec := gob.NewDecoder(reader)
	if err := dec.Decode(&parameters); err != nil {
		return err
	}

	context.SetParameters(parameters.N, parameters.Modulus)
	context.ValidateParameters()

	return nil
}


// IsValidated returns true if the context has been validated (for NTT), else false.
func (context *Context) IsValidated() bool {
	return context.validated
}

// GetBRedParams returns the Barret reduction parameters of the context.
func (context *Context) GetBredParams() [][]uint32 {
	return context.bredParams
}

// GetMredParams returns the Montgomery reduction parameters of the context.
func (context *Context) GetMredParams() []uint32 {
	return context.mredParams
}

// GetPsi returns the primitive root used to compute the NTT parameters of the context.
func (context *Context) GetPsi() []uint32 {
	return context.psiMont
}

// GetPsi returns the primitive root used to compute the InvNTT parameters of the context.
func (context *Context) GetPsiInv() []uint32 {
	return context.psiInvMont
}

// GetNttPsi returns the NTT parameters of the context.
func (context *Context) GetNttPsi() [][]uint32 {
	return context.nttPsi
}

//GetNttPsiInv returns the InvNTT parameters of the context.
func (context *Context) GetNttPsiInv() [][]uint32 {
	return context.nttPsiInv
}

// GetNttNInv returns 1/N mod each moduli.
func (context *Context) GetNttNInv() []uint32 {
	return context.nttNInv
}

// NewPoly create a new polynomial with all coefficients set to 0.
func (context *Context) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint32, len(context.Modulus))
	for i := 0; i < len(context.Modulus); i++ {
		p.Coeffs[i] = make([]uint32, context.N)
	}

	return p
}

// NewUniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1]
func (context *Context) NewUniformPoly() (Pol *Poly) {

	var randomBytes []byte
	var randomUint, mask uint32

	Pol = context.NewPoly()

	n := context.N
	if n < 4 {
		n = 4
	}

	randomBytes = make([]byte, n)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	for j, qi := range context.Modulus {

		// Starts by computing the mask
		mask = (1 << uint32(bits.Len32(qi))) - 1

		// Iterates for each modulus over each coefficient
		for i := uint32(0); i < context.N; i++ {

			// Samples an integer between [0, qi-1]
			for {

				// Replenishes the pool if it runs empty
				if len(randomBytes) < 4 {
					randomBytes = make([]byte, n)
					if _, err := rand.Read(randomBytes); err != nil {
						panic("crypto rand error")
					}
				}

				// Reads bytes from the pool
				randomUint = binary.BigEndian.Uint32(randomBytes[:4]) & mask
				randomBytes = randomBytes[4:] // Discard the used bytes

				// If the integer is between [0, qi-1], breaks the loop
				if randomUint < qi {
					break
				}
			}

			Pol.Coeffs[j][i] = randomUint
		}
	}

	return
}
