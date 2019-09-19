package newhope

import (
	"crypto/rand"
	"encoding/binary"
	"errors"
	"math/bits"
)

//==============================
//===== POLYNOMIAL CONTEXT =====
//==============================

// Context is a structure keeping all the variable required to operate on a polynomial represented in this context.
// This include its modulus, modular reduction and ntt transformation.
type Context struct {

	// Polynomial nb.Coefficients
	n uint32

	// Modulus
	q uint32

	// 2^bit_length(Qi) - 1
	mask uint32

	// Determines if NTT can be used with the current context.
	validated bool

	// Fast reduction parameters
	bredparam uint64
	mredparam uint32

	//NTT Parameters
	psiMont    uint32 //2nth primitive root in montgomery form
	psiInvMont uint32 //2nth inverse primitive root in montgomery form

	nttPsi    []uint32 //powers of the inverse of the 2nth primitive root in montgomery form (in bitreversed order)
	nttPsiInv []uint32 //powers of the inverse of the 2nth primitive root in montgomery form (in bitreversed order)
	nttNInv   uint32   //[N^-1] mod q in montgomery form
}

// NewContext generates a new empty context.
func NewContext(N uint32, Q uint32) (context *Context, err error) {
	context = new(Context)

	if err = context.setparameters(N, Q) ; err != nil {
		return nil, err
	}

	if err = context.validateparameters() ; err != nil {
		return context, err
	}

	return
}


// SetParameters initialize the parameters of an empty context with N and the provided moduli.
// Only checks that N is a power of 2 and computes all the variable that aren't used for the NTT.
func (context *Context) setparameters(N uint32, Q uint32) error {

	// Checks if N is a power of 2
	if (N&(N-1)) != 0 && N != 0 {
		return errors.New("invalid ring degree (must be a power of 2)")
	}

	context.validated = false

	context.n = N

	context.q = Q

	context.mask = (1 << uint32(bits.Len32(Q))) - 1

	//Computes the fast modular reduction parameters for the Context
	context.bredparam = bredparam(Q)

	// If qi is not a power of 2, we can compute the mredparams (else it should not
	// because it will return an error and there is no valid montgomery form mod a power of 2)
	if (Q&(Q-1)) != 0 && Q != 0 {
		context.mredparam = mredparam(Q)
	}

	return nil
}

// ValidateParameters checks that N has beed correctly initialized, and checks that the modulus is a prime congruent to 1 mod 2N (i.e. allowing NTT).
// Then it computes the variables required for the NTT. ValidateParameters purpose is to validate that the moduli allow the NTT and compute the
// NTT parameters.
func (context *Context) validateparameters() error {

	if context.validated {
		return nil
	}

	if context.n == 0 || context.q == 0 {
		return errors.New("error : invalid context parameters (missing)")
	}

	// CHECKS IF VALIDE NTT
	// Checks if each qi is Prime and if qi = 1 mod 2n

	if isprime(context.q) == false || context.q&((context.n<<1)-1) != 1 {
		context.validated = false
		return errors.New("warning : provided q does not allow NTT")
	}

	context.nttPsi = make([]uint32, context.n)
	context.nttPsiInv = make([]uint32, context.n)

	bitLenofN := uint32(bits.Len32(context.n) - 1)

	//2.1 Computes N^(-1) mod Q in Montgomery form
	context.nttNInv = mform(modexp(context.n, context.q-2, context.q), context.q, context.bredparam)

	//2.2 Computes Psi and PsiInv in Montgomery form
	context.nttPsi = make([]uint32, context.n)
	context.nttPsiInv = make([]uint32, context.n)

	//Finds a 2nth primitive Root
	g := primitiveroot(context.q)

	_2n := uint32(context.n << 1)

	power := (context.q - 1) / _2n
	powerInv := (context.q - 1) - power

	//Computes Psi and PsiInv in Montgomery Form
	PsiMont := mform(modexp(g, power, context.q), context.q, context.bredparam)
	PsiInvMont := mform(modexp(g, powerInv, context.q), context.q, context.bredparam)

	context.psiMont = PsiMont
	context.psiInvMont = PsiInvMont

	context.nttPsi[0] = mform(1, context.q, context.bredparam)
	context.nttPsiInv[0] = mform(1, context.q, context.bredparam)

	// Computes nttPsi[j] = nttPsi[j-1]*Psi and nttPsiInv[j] = nttPsiInv[j-1]*PsiInv
	for j := uint32(1); j < context.n; j++ {

		indexReversePrev := bitreverse32(j-1, bitLenofN)
		indexReverseNext := bitreverse32(j, bitLenofN)

		context.nttPsi[indexReverseNext] = mred(context.nttPsi[indexReversePrev], PsiMont, context.q, context.mredparam)
		context.nttPsiInv[indexReverseNext] = mred(context.nttPsiInv[indexReversePrev], PsiInvMont, context.q, context.mredparam)
	}

	context.validated = true

	return nil
}

// N returns the ring degree
func (context *Context) N() uint32 {
	return context.n
}

// Modulus returns the ring modulus
func (context *Context) Modulus() uint32 {
	return context.q
}

// IsValidated returns true if the context has been validated (for NTT), else false.
func (context *Context) IsValidated() bool {
	return context.validated
}

// Getbredparam returns the Barret reduction parameters of the context.
func (context *Context) Getbredparam() uint64 {
	return context.bredparam
}

// Getmredparams returns the Montgomery reduction parameters of the context.
func (context *Context) Getmredparams() uint32 {
	return context.mredparam
}

// GetPsi returns the primitive root used to compute the NTT parameters of the context.
func (context *Context) GetPsi() uint32 {
	return context.psiMont
}

// GetPsi returns the primitive root used to compute the InvNTT parameters of the context.
func (context *Context) GetPsiInv() uint32 {
	return context.psiInvMont
}

// GetNttPsi returns the NTT parameters of the context.
func (context *Context) GetNttPsi() []uint32 {
	return context.nttPsi
}

//GetNttPsiInv returns the InvNTT parameters of the context.
func (context *Context) GetNttPsiInv() []uint32 {
	return context.nttPsiInv
}

// GetNttNInv returns 1/N mod each moduli.
func (context *Context) GetNttNInv() uint32 {
	return context.nttNInv
}

// NewPoly create a new polynomial with all coefficients set to 0.
func (context *Context) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([]uint32, context.n)

	return p
}

// NewUniformPoly generates a new polynomial with coefficients following a uniform distribution over [0, q-1]
func (context *Context) NewUniformPoly() (Pol *Poly) {

	var randomBytes []byte
	var randomUint uint32

	Pol = context.NewPoly()

	n := context.n
	if n < 4 {
		n = 4
	}

	randomBytes = make([]byte, n)
	if _, err := rand.Read(randomBytes); err != nil {
		panic("crypto rand error")
	}

	// Iterates for each q over each coefficient
	for i := uint32(0); i < context.n; i++ {

		// Samples an integer between [0, q-1]
		for {

			// Replenishes the pool if it runs empty
			if len(randomBytes) < 4 {
				randomBytes = make([]byte, n)
				if _, err := rand.Read(randomBytes); err != nil {
					panic("crypto rand error")
				}
			}

			// Reads bytes from the pool
			randomUint = binary.BigEndian.Uint32(randomBytes[:4]) & context.mask
			randomBytes = randomBytes[4:] // Discard the used bytes

			// If the integer is between [0, q-1], breaks the loop
			if randomUint < context.q {
				break
			}
		}

		Pol.Coeffs[i] = randomUint
	}

	return
}
