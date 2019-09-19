package newhope

import (
	"fmt"
	"math/rand"
	"testing"
	"time"
)

func Benchmark_newhope(b *testing.B) {

	rand.Seed(time.Now().UnixNano())

	var N, Q uint32

	N = 1024
	Q = 12289

	//sigma := 3.19

	context := NewContext()
	context.SetParameters(N, Q)
	context.ValidateParameters()

	benchmark_BRed(context, b)
	benchmark_BRedAdd(context, b)
	benchmark_MRed(context, b)
	benchmark_KRed(context, b)
	//benchmark_KRed2x(context, b)

}

func benchmark_BRed(context *Context, b *testing.B) {

	x := rand.Uint32() % context.Modulus
	y := rand.Uint32() % context.Modulus

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%d/BRed", context.Modulus), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = BRed(x, y, context.Modulus, context.bredParams)
		}
	})
}

func benchmark_BRedAdd(context *Context, b *testing.B) {

	x := rand.Uint32()

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%dbit/BRedAdd", context.Modulus), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BRedAdd(x, context.Modulus, context.bredParams)
		}
	})
}

func benchmark_MRed(context *Context, b *testing.B) {

	x := rand.Uint32() % context.Modulus
	y := rand.Uint32() % context.Modulus

	y = MForm(y, context.Modulus, context.bredParams)

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%d/MRed", context.Modulus), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = MRed(x, y, context.Modulus, context.mredParams)
		}
	})
}

func benchmark_KRed(context *Context, b *testing.B) {

	x := int64(rand.Uint32() % context.Modulus)
	y := int64(rand.Uint32() % context.Modulus)

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%d/KRed", context.Modulus), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = KRed(x*y)
		}
	})
}


