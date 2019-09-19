package newhope

import (
	"fmt"
	"math/rand"
	"testing"
	"time"
	"log"
)

func Test_newhope(t *testing.T) {

	var err error

	rand.Seed(time.Now().UnixNano())

	var N, Q uint32

	N = 1024
	Q = 12289

	//sigma := 3.19

	var context *Context
	if context, err = NewContext(N, Q) ; err != nil {
		log.Fatal(err)
	}

	test_BRed(context, t)
	test_BRedAdd(context, t)
	test_MRed(context, t)

	test_MulPoly(context, t)
	test_MulPoly_Montgomery(context, t)

}

func test_BRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/bred", context.n, context.q), func(t *testing.T) {

		var x, y uint32
		var z uint64

		for i := 0; i < 0x10000; i++ {
			x = rand.Uint32() % context.q
			y = rand.Uint32() % context.q
			z = (uint64(x) * uint64(y)) % uint64(context.q)

			test := bred(x, y, context.q, context.bredparam)

			if test != uint32(z) {
				t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, z)
				break
			}
		}
	})
}

func test_BRedAdd(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/bredadd", context.n, context.q), func(t *testing.T) {

		var x, y uint32
		var z uint64

		for i := 0; i < 0x10000; i++ {
			x = rand.Uint32() % context.q
			y = rand.Uint32() % context.q
			z = (uint64(x) + uint64(y)) % uint64(context.q)

			test := bredadd(x+y, context.q, context.bredparam)

			if test != uint32(z) {
				t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, z)
				break
			}
		}
	})
}

func test_MRed(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/mred", context.n, context.q), func(t *testing.T) {

		var x, y uint32
		var z uint64

		for i := 0; i < 0x10000; i++ {
			x = rand.Uint32() % context.q
			y = rand.Uint32() % context.q
			z = (uint64(x) * uint64(y)) % uint64(context.q)

			test := mred(x, mform(y, context.q, context.bredparam), context.q, context.mredparam)

			if test != uint32(z) {
				t.Errorf("error : 128bit barrett multiplication, x = %v, y=%v, have = %v, want =%v", x, y, test, z)
				break
			}
		}
	})
}

func test_MulPoly(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/MulPoly", context.n, context.q), func(t *testing.T) {

		p1 := context.NewUniformPoly()
		p2 := context.NewUniformPoly()

		p3Test := context.NewPoly()
		p3Want := context.NewPoly()

		context.MulPolyNaive(p1, p2, p3Want)
		context.MulPoly(p1, p2, p3Test)

		for i := uint32(0); i < context.n; i++ {
			if p3Want.Coeffs[i] != p3Test.Coeffs[i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[i], p3Test.Coeffs[i])
				break
			}
		}
	})
}

func test_MulPoly_Montgomery(context *Context, t *testing.T) {

	t.Run(fmt.Sprintf("N=%d/Q=%d/MulPoly_Montgomery", context.n, context.q), func(t *testing.T) {
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

		for i := uint32(0); i < context.n; i++ {
			if p3Want.Coeffs[i] != p3Test.Coeffs[i] {
				t.Errorf("ERROR MUL COEFF %v, want %v - has %v", i, p3Want.Coeffs[i], p3Test.Coeffs[i])
				break
			}
		}
	})
}
