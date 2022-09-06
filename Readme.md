# A rust implementation of a Quadratic Sieve 
(it is not a Quadratic Sieve but nearly, it's more like a Dixon's factorization method)


## Exemple 

```
./quadratic_sieve (master) » time cargo run --release
    Finished release [optimized] target(s) in 0.01s
     Running `target/release/quadratic_sieve`
[2022-09-06T17:29:22.439Z INFO  quadratic_sieve] Let's try to find factors for N=18070820886874048059516561644059055662781
...
[2022-09-06T18:27:12.200Z WARN  quadratic_sieve] [733, 210109, 117335452064317970808827705689973] product: 18070820886874048059516561644059055662781 (n = 18070820886874048059516561644059055662781)
cargo run --release  41621.20s user 2.26s system 1199% cpu 57:49.81 total
```

## References: 
- [Wikipedia (Quadratic Sieve)](https://en.wikipedia.org/wiki/Quadratic_sieve)
- [Wikipedia (Dixon's factorization methode)](https://en.wikipedia.org/wiki/Dixon%27s_factorization_method)
- [The Quadratic Sieve Factoring Algorithm](https://www.cs.virginia.edu/crab/QFS_Simple.pdf)
