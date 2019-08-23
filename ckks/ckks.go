package ckks

import (
	"errors"
	"github.com/lca1/lattigo/ring"
	"math"
	"math/bits"
)

type CkksContext struct {

	// Context parameters
	logN     uint64
	logQ     uint64
	logScale uint64
	n        uint64
	slots    uint64

	// Uperbound in bits of the modulie
	maxBit uint64

	// Number of avaliable levels
	levels uint64

	// Modulie chain
	modulie []uint64

	// Contexts chain
	contextLevel []*ring.Context

	// Keys' context
	keyscontext *ring.Context

	// Pre-computed values for the rescaling
	scalechain   []uint64
	rescalParams [][]uint64

	// Sampling variance
	sigma float64

	// Samplers
	gaussianSampler *ring.KYSampler
	ternarySampler  *ring.TernarySampler

	// Galoi generator for the rotations, encoding and decoding params
	gen    uint64
	genInv uint64

	// Rotation params
	galElRotRow      uint64
	galElRotColLeft  []uint64
	galElRotColRight []uint64

	// Encoding and Decoding params
	indexMatrix []uint64
	gap         uint64
	roots       []complex128
	inv_roots   []complex128

	// Checksum of [n, [modulies]]
	checksum []byte
}

// NewCkksContext generates a new ckkscontext, given the parameters logN, logQ, logScale (in base 2),
// levels and sigma.
func NewCkksContext(logN uint64, moduliechain []uint64, logScale uint64, sigma float64) (*CkksContext, error) {

	var err error

	ckkscontext := new(CkksContext)

	ckkscontext.logN = logN
	ckkscontext.n = 1 << logN
	ckkscontext.slots = 1 << (logN - 1)
	ckkscontext.logScale = logScale
	ckkscontext.levels = uint64(len(moduliechain))

	// ========== START < PRIMES GENERATION > START ===============
	ckkscontext.scalechain = make([]uint64, len(moduliechain))

	primesbitlen := make(map[uint64]uint64)

	for i := range moduliechain {

		primesbitlen[moduliechain[i]] += 1

		ckkscontext.scalechain[i] = moduliechain[i]

		if moduliechain[i] > 60 {
			return nil, errors.New("error : provided modulie must be smaller than 60")
		}
	}

	primes := make(map[uint64][]uint64)

	for key, value := range primesbitlen {
		primes[key], _ = GenerateCKKSPrimes(key, logN, value)
	}

	ckkscontext.modulie = make([]uint64, len(moduliechain))

	for i := range moduliechain {
		ckkscontext.modulie[i] = primes[moduliechain[i]][0]
		primes[moduliechain[i]] = primes[moduliechain[i]][1:]

		if uint64(bits.Len64(ckkscontext.modulie[i])) > ckkscontext.maxBit {
			ckkscontext.maxBit = uint64(bits.Len64(ckkscontext.modulie[i]))
		}
	}

	// ========== END < PRIMES GENERATION > END ===============

	// ========== START < CONTEXTS CHAIN > START ===============
	ckkscontext.contextLevel = make([]*ring.Context, ckkscontext.levels)

	ckkscontext.contextLevel[0] = ring.NewContext()

	if err = ckkscontext.contextLevel[0].SetParameters(1<<logN, ckkscontext.modulie[:1]); err != nil {
		return nil, err
	}

	if err = ckkscontext.contextLevel[0].ValidateParameters(); err != nil {
		return nil, err
	}

	for i := uint64(1); i < ckkscontext.levels; i++ {

		ckkscontext.contextLevel[i] = ring.NewContext()

		if err = ckkscontext.contextLevel[i].SetParameters(1<<logN, ckkscontext.modulie[i:i+1]); err != nil {
			return nil, err
		}

		if err = ckkscontext.contextLevel[i].ValidateParameters(); err != nil {
			return nil, err
		}

		// Instead of recomputing and storing redundant context, subsequent contexts are a merge of previous contexts.
		if err = ckkscontext.contextLevel[i].Merge(ckkscontext.contextLevel[i-1], ckkscontext.contextLevel[i]); err != nil {
			return nil, err
		}
	}
	// ========== END < CONTEXTS CHAIN > END ===============

	// Context used for the generation of the keys
	ckkscontext.keyscontext = ckkscontext.contextLevel[ckkscontext.levels-1]
	ckkscontext.logQ = uint64(ckkscontext.contextLevel[ckkscontext.levels-1].ModulusBigint.Value.BitLen())

	// ========== START < RESCALE PRE-COMPUATION PARAMETERS > START ===============
	var Qi, Ql uint64

	ckkscontext.rescalParams = make([][]uint64, ckkscontext.levels-1)

	for j := ckkscontext.levels - 1; j > 0; j-- {

		ckkscontext.rescalParams[j-1] = make([]uint64, j)

		Ql = ckkscontext.modulie[j]

		bredParams := ckkscontext.contextLevel[j-1].GetBredParams()

		for i := uint64(0); i < j; i++ {

			Qi = ckkscontext.modulie[i]

			ckkscontext.rescalParams[j-1][i] = ring.MForm(modexp(Ql, Qi-2, Qi), Qi, bredParams[i])
		}
	}
	// ========== END < RESCALE PRE-COMPUATION PARAMETERS > END ===============

	// default variance
	ckkscontext.sigma = sigma

	ckkscontext.gaussianSampler = ckkscontext.keyscontext.NewKYSampler(sigma, int(6*sigma))
	ckkscontext.ternarySampler = ckkscontext.keyscontext.NewTernarySampler()

	// ========== START < ROTATION ELEMENTS > START ===============
	var m, mask uint64

	m = ckkscontext.n << 1

	mask = m - 1

	ckkscontext.gen = 5 // Any integer equal to 1 mod 4 and comprime to 2N will do fine
	ckkscontext.genInv = modexp(ckkscontext.gen, mask, m)

	ckkscontext.galElRotColLeft = make([]uint64, ckkscontext.slots)
	ckkscontext.galElRotColRight = make([]uint64, ckkscontext.slots)

	ckkscontext.galElRotColRight[0] = 1
	ckkscontext.galElRotColLeft[0] = 1

	for i := uint64(1); i < ckkscontext.slots; i++ {
		ckkscontext.galElRotColLeft[i] = (ckkscontext.galElRotColLeft[i-1] * ckkscontext.gen) & mask
		ckkscontext.galElRotColRight[i] = (ckkscontext.galElRotColRight[i-1] * ckkscontext.genInv) & mask
	}

	ckkscontext.galElRotRow = mask
	// ============ END < ROTATION ELEMENTS > END ================

	// ========== START < ENCODER PARAMETERS > START ================
	var pos, index1, index2 uint64

	ckkscontext.gap = 1 //gap-1 is the gap between each slot, here 1 means no gap.

	ckkscontext.indexMatrix = make([]uint64, ckkscontext.n)

	pos = 1

	for i := uint64(0); i < ckkscontext.slots; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		ckkscontext.indexMatrix[i] = bitReverse64(index1, ckkscontext.logN)
		ckkscontext.indexMatrix[i|ckkscontext.slots] = bitReverse64(index2, ckkscontext.logN)

		pos *= ckkscontext.gen
		pos &= mask
	}

	ckkscontext.roots = make([]complex128, ckkscontext.n)
	ckkscontext.inv_roots = make([]complex128, ckkscontext.n)

	angle := 6.283185307179586 / float64(m)
	psi := complex(math.Cos(angle), math.Sin(angle))
	psiInv := complex(1, 0) / psi

	ckkscontext.roots[0] = 1
	ckkscontext.inv_roots[0] = 1

	for j := uint64(1); j < ckkscontext.n; j++ {

		indexReversePrev := bitReverse64(j-1, ckkscontext.logN)
		indexReverseNext := bitReverse64(j, ckkscontext.logN)

		ckkscontext.roots[indexReverseNext] = ckkscontext.roots[indexReversePrev] * psi
		ckkscontext.inv_roots[indexReverseNext] = ckkscontext.inv_roots[indexReversePrev] * psiInv
	}
	// ========== END < ENCODER PARAMETERS > END ================

	// =========== START < CHECKSUM > =============
	toHash := make([]uint64, 1+ckkscontext.levels)
	toHash[0] = logN
	for i := 1; i < len(toHash); i++ {
		toHash[i] = ckkscontext.modulie[i-1]
	}

	if ckkscontext.checksum, err = Hash(toHash); err != nil {
		return nil, err
	}
	// ========== END < CHECKSUM > END

	return ckkscontext, nil

}

func (ckkscontext *CkksContext) LogN() uint64 {
	return ckkscontext.logN
}

func (ckkscontext *CkksContext) LogQ() uint64 {
	return ckkscontext.logQ
}

func (ckkscontext *CkksContext) Modulie() []uint64 {
	return ckkscontext.modulie
}

func (ckkscontext *CkksContext) Levels() uint64 {
	return ckkscontext.levels
}

func (ckkscontext *CkksContext) LogScale() uint64 {
	return ckkscontext.logScale
}

func (ckkscontext *CkksContext) ContextKeys() *ring.Context {
	return ckkscontext.keyscontext
}

func (ckkscontext *CkksContext) Slots() uint64 {
	return (1 << (ckkscontext.logN - 1))
}
