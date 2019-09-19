package newhope

import (
	"fmt"
	"math/rand"
	"testing"
	"time"
)

func Test_newhope(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	var N, Q uint32

	N = 4
	Q = 12289

	//sigma := 3.19

	context := NewContext()
	context.SetParameters(N, Q)
	context.ValidateParameters()

	test_BRed(context, t)
	test_BRedAdd(context, t)
	test_MRed(context, t)

	test_MulPoly(context, t)
	test_MulPoly_Montgomery(context, t)

}

func test_BRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/BRed", context.N, context.Modulus), func(t *testing.T) {

		var x, y uint32
		var z uint64

		for i := 0; i < 0x10000; i++ {
			x = rand.Uint32() % context.Modulus
			y = rand.Uint32() % context.Modulus
			z = (uint64(x) * uint64(y)) % uint64(context.Modulus)

			test := BRed(x, y, context.Modulus, context.bredParams)

			if test != uint32(z) {
				t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, z)
				break
			}
		}
	})
}

func test_BRedAdd(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/BRedAdd", context.N, context.Modulus), func(t *testing.T) {

		var x, y uint32
		var z uint64

		for i := 0; i < 0x10000; i++ {
			x = rand.Uint32() % context.Modulus
			y = rand.Uint32() % context.Modulus
			z = (uint64(x) + uint64(y)) % uint64(context.Modulus)

			test := BRedAdd(x+y, context.Modulus, context.bredParams)

			if test != uint32(z) {
				t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, z)
				break
			}
		}
	})
}

func test_MRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/MRed", context.N, context.Modulus), func(t *testing.T) {

		var x, y uint32
		var z uint64

		for i := 0; i < 0x10000; i++ {
			x = rand.Uint32() % context.Modulus
			y = rand.Uint32() % context.Modulus
			z = (uint64(x) * uint64(y)) % uint64(context.Modulus)

			test := MRed(x, MForm(y, context.Modulus, context.bredParams), context.Modulus, context.mredParams)

			if test != uint32(z) {
				t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, z)
				break
			}
		}
	})
}

func test_MulPoly(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/MulPoly", context.N, context.Modulus), func(t *testing.T) {

		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()

		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MulPolyNaive(p1, p2, p3Want)
		context.MulPoly(p1, p2, p3Test)

		for i := uint32(0); i < context.N; i++ {
			if p3Want.Coeffs[i] != p3Test.Coeffs[i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[i], p3Test.Coeffs[i])
				break
			}
		}
	})
}

func test_MulPoly_Montgomery(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/MulPoly_Montgomery", context.N, context.Modulus), func(t *testing.T) {
		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()

		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MForm(p1, p1)
		context.MForm(p2, p2)

		context.MulPolyNaiveMontgomery(p1, p2, p3Want)
		context.MulPolyMontgomery(p1, p2, p3Test)

		context.InvMForm(p3Test, p3Test)
		context.InvMForm(p3Want, p3Want)

		for i := uint32(0); i < context.N; i++ {
			if p3Want.Coeffs[i] != p3Test.Coeffs[i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[i], p3Test.Coeffs[i])
				break
			}
		}
	})
}
