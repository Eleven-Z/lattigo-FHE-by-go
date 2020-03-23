//Shares-to-encryption is the protocol described in section IV.E.2 of the paper.
//Unlike the encryption-to-shares, no distinction is made between nodes: there is no master and slaves.
//Nodes first need to sample the CRS, which will be the degree-1 term of the final ciphertext.
//With the CRS, they produce a re-encryption share, which is an ordinary CKS share (with sk_input = 0)
//to which they sum delta*M_i (M_i being their initial additive share, in Z_t^n).
//The degree-0 term of the ciphertext is then the sum of the re-encryption shares.

package dbfv

import (
	"github.com/ldsec/lattigo/bfv"
	"github.com/ldsec/lattigo/ring"
)

//S2EProtocol contains all the parameters needed to perform the various steps of the protocol.
type S2EProtocol struct {
	cks     *CKSProtocol //CKSProtocol is not embedded to have control over the exposed methods
	encoder bfv.Encoder

	//Memory pools
	plain  *bfv.Plaintext
	cipher *bfv.Ciphertext
}

//S2EReencryptionShare represents the re-encryption share. They need to be all collected and aggregated.
//It is an element of R_q.
type S2EReencryptionShare struct {
	CKSShare
}

//NewS2EProtocol allocates a protocol struct.
func NewS2EProtocol(params *bfv.Parameters, sigmaSmudging float64) *S2EProtocol {
	cks := NewCKSProtocol(params, sigmaSmudging)

	return &S2EProtocol{cks,
		bfv.NewEncoder(params),
		bfv.NewPlaintext(params),
		bfv.NewCiphertext(params, 1)}
}

//AllocateShare allocates a re-encryption share.
func (s2e *S2EProtocol) AllocateShare() *S2EReencryptionShare {
	return &S2EReencryptionShare{s2e.cks.AllocateShare()}
}

//GenShare generates the re-encryption share, given the output secret key, the additive share, and the crs.
func (s2e *S2EProtocol) GenShare(sk *bfv.SecretKey, crs *ring.Poly, addShare *AdditiveShare, shareOut *S2EReencryptionShare) {
	//First step is to run the CKS protocol, on a crafted ciphertext, with s_in = 0.
	//First, craft the ciphertext.
	s2e.encoder.EncodeUint(addShare.elem.GetCoefficients()[0], s2e.plain) //Store delta*M_i, which will be ct[0]
	s2e.cipher.SetValue([]*ring.Poly{s2e.plain.Value()[0], crs})          //Build the ciphertext
	//Then, run the CKS protocol
	s2e.cks.GenShare(s2e.cks.context.contextQ.NewPoly(), sk.Get(), s2e.cipher, shareOut.CKSShare)

	//We add delta*M_i to the re-encryption share
	s2e.cks.context.contextQ.Add(shareOut.Poly, s2e.plain.Value()[0], shareOut.Poly)

	return
}

//AggregateShares pretty much describes itself. It is safe to have shareOut coincide with share1 or share2.
func (s2e *S2EProtocol) AggregateShares(share1, share2, shareOut *S2EReencryptionShare) {
	s2e.cks.context.contextQ.Add(share1.Poly, share2.Poly, shareOut.Poly)
}

//Reencrypt takes the aggregate of all re-encryption share and the crs to build a (fresh) encryption
func (s2e *S2EProtocol) Reencrypt(shareAgg *S2EReencryptionShare, crs *ring.Poly) (ct *bfv.Ciphertext) {
	ct = bfv.NewCiphertext(s2e.cks.context.params, 1)
	ct.SetValue([]*ring.Poly{shareAgg.Poly, crs})

	return
}

/******** Operations on additive shares********/

func (s2e *S2EProtocol) GenRandomAddShare() *AdditiveShare {
	poly := s2e.cks.context.contextT.NewUniformPoly()
	return &AdditiveShare{poly}
}
