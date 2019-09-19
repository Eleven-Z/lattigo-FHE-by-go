package newhope

import (
	"fmt"
	"math/rand"
	"testing"
	"time"
	"log"
)

func Benchmark_newhope(b *testing.B) {

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

	benchmark_BRed(context, b)
	benchmark_BRedAdd(context, b)
	benchmark_MRed(context, b)
	benchmark_KRed(context, b)
	//benchmark_KRed2x(context, b)

}

func benchmark_BRed(context *Context, b *testing.B) {

	x := rand.Uint32() % context.q
	y := rand.Uint32() % context.q

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%d/bred", context.q), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = bred(x, y, context.q, context.bredparam)
		}
	})
}

func benchmark_BRedAdd(context *Context, b *testing.B) {

	x := rand.Uint32()

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%dbit/bredadd", context.q), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			bredadd(x, context.q, context.bredparam)
		}
	})
}

func benchmark_MRed(context *Context, b *testing.B) {

	x := rand.Uint32() % context.q
	y := rand.Uint32() % context.q

	y = mform(y, context.q, context.bredparam)

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%d/MRed", context.q), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = mred(x, y, context.q, context.mredparam)
		}
	})
}

func benchmark_KRed(context *Context, b *testing.B) {

	x := int64(rand.Uint32() % context.q)
	y := int64(rand.Uint32() % context.q)

	b.ResetTimer()

	b.Run(fmt.Sprintf("Q=%d/KRed", context.q), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = kred(x*y)
		}
	})
}


