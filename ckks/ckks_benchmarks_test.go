package ckks

import (
	"fmt"
	"testing"
)

type benchParams struct {
	logN     uint64
	moduli   []uint64
	logScale uint64
	sigma    float64
	bdc      uint64
}

func BenchmarkCKKSScheme(b *testing.B) {

	var err error
	var ckkscontext *CkksContext
	var kgen *KeyGenerator
	var sk *SecretKey
	var pk *PublicKey
	var rlk *EvaluationKey
	var rotkey *RotationKey
	var encryptor *Encryptor
	var decryptor *Decryptor
	var evaluator *Evaluator
	var plaintext *Plaintext
	var ciphertext *Ciphertext

	params := []benchParams{
		//{logN: 10, logQ: 30, levels:  1, logScale : 20, sigma: 3.19, bdc : 60},
		//{logN: 11, logQ: 60, levels:  1, logScale : 40, sigma: 3.19, bdc : 60},
		{logN: 12, moduli: []uint64{40, 40}, logScale: 40, sigma: 3.2, bdc: 60},
		{logN: 13, moduli: []uint64{40, 40, 40, 40}, logScale: 40, sigma: 3.2, bdc: 60},
		{logN: 14, moduli: []uint64{40, 40, 40, 40, 40, 40, 40, 40}, logScale: 40, sigma: 3.2, bdc: 60},
		{logN: 15, moduli: []uint64{40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40}, logScale: 40, sigma: 3.2, bdc: 60},
	}

	var logN, logQ, logScale, levels, bdc uint64
	var sigma float64

	for _, param := range params {

		logN = param.logN
		logScale = param.logScale
		sigma = param.sigma
		bdc = param.bdc

		if ckkscontext, err = NewCkksContext(logN, param.moduli, logScale, sigma); err != nil {
			b.Error(err)
		}

		kgen = ckkscontext.NewKeyGenerator()

		if sk, pk, err = kgen.NewKeyPair(); err != nil {
			b.Error(err)
		}

		if rlk, err = kgen.NewRelinKey(sk, bdc); err != nil {
			b.Error(err)
		}

		if rotkey, err = kgen.NewRotationKeysPow2(sk, bdc, true); err != nil {
			b.Error(err)
		}

		if encryptor, err = ckkscontext.NewEncryptor(pk); err != nil {
			b.Error(err)
		}

		if decryptor, err = ckkscontext.NewDecryptor(sk); err != nil {
			b.Error(err)
		}

		evaluator = ckkscontext.NewEvaluator()

		_ = rotkey

		ciphertext = ckkscontext.NewCiphertext(1, levels-1, logScale)

		var values []complex128
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/Levels=%d/decomp=%d/sigma=%.2f/Encode", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			slots := 1 << (logN - 2)
			values = make([]complex128, slots)
			for i := 0; i < slots; i++ {
				values[i] = complex(randomFloat(0.1, 1), 0)
			}

			plaintext = ckkscontext.NewPlaintext(ckkscontext.levels-1, logScale)

			for i := 0; i < b.N; i++ {
				if err = plaintext.EncodeComplex(values); err != nil {
					b.Error(err)
				}
			}
		})

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Decode", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if values = plaintext.DecodeComplex(); err != nil {
					b.Error(err)
				}
			}
		})

		// Key Pair Generation
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/KeyPairGen", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				sk, pk, err = kgen.NewKeyPair()
				if err != nil {
					b.Error(err)
				}
			}
		})

		// SwitchKeyGen
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/SwitchKeyGen", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if rlk, err = kgen.NewRelinKey(sk, bdc); err != nil {
					b.Error(err)
				}
			}
		})

		// Encrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Encrypt", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = encryptor.Encrypt(plaintext, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Decrypt
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Decrypt", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = decryptor.Decrypt(ciphertext, plaintext); err != nil {
					b.Error(err)
				}
			}
		})

		// Add
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Add", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Add(ciphertext, ciphertext, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Add Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/AddScalar", logN, logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.AddConst(ciphertext, complex(3.1415, -1.4142), ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Mult Scalar
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/MultScalar", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MultConst(ciphertext, complex(3.1415, -1.4142), ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul
		receiver := ckkscontext.NewRandomCiphertext(2, ckkscontext.levels-1, logScale)
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Multiply", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext, ciphertext, nil, receiver); err != nil {
					b.Error(err)
				}
			}
		})

		// Mul Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/MulRelin", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.MulRelin(ciphertext, ciphertext, rlk, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Square
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Square", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Square(ciphertext, nil, receiver); err != nil {
					b.Error(err)
				}
			}
		})

		// Square Relin
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/SquareRelin", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Square(ciphertext, rlk, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Relin

		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Relin", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Relinearize(receiver, rlk, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Rescale
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Rescale", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				rescale(evaluator, ciphertext.Value()[0], receiver.Value()[0])
				rescale(evaluator, ciphertext.Value()[1], receiver.Value()[1])
			}
		})

		// Conjugate / Rotate
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/Conjugate", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.Conjugate(ciphertext, rotkey, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})

		// Rotate Cols
		b.Run(fmt.Sprintf("logN=%d/logQ=%d/levels=%d/decomp=%d/sigma=%.2f/RotateCols", logN, ckkscontext.logQ, ckkscontext.levels, bdc, sigma), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if err = evaluator.RotateColumns(ciphertext, 1, rotkey, ciphertext); err != nil {
					b.Error(err)
				}
			}
		})
	}
}
