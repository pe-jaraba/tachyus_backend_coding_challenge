let rec _subtract (u: array<float>, v: array<float>, i: int, accum: array<float>) = 
    if i = u.Length then
        accum
    else
        Array.set accum i (u.[i] - v.[i])
        _subtract(u, v, i + 1, accum)

let subtract (u: array<float>, v: array<float>) =
    _subtract(u, v, 0, Array.zeroCreate<float> u.Length)

let rec _sum (u: array<float>, v: array<float>, i: int, accum: array<float>) = 
    if i = u.Length then
        accum
    else
        Array.set accum i (u.[i] + v.[i])
        _sum(u, v, i + 1, accum)

let sum (u: array<float>, v: array<float>) =
    _sum(u, v, 0, Array.zeroCreate<float> u.Length)

let rec _scale (x: float, v: array<float>, i: int, accum: array<float>) =
    if i = v.Length then
        accum
    else
        Array.set accum i (x*v.[i])
        _scale(x, v, i + 1, accum)

let scale (x: float, v: array<float>) =
    _scale(x, v, 0, Array.zeroCreate<float> v.Length)

let rec _multiplyAt (V: array<array<float>>, index: int, i: int, accum: float) =
    if i = V.Length then
        accum
    else
        if i = 0 then
            _multiplyAt(V, index, i + 1, V.[i].[index])
        else
            _multiplyAt(V, index, i + 1, accum * V.[i].[index])

let multiplyAt (V: array<array<float>>, index: int) =
    _multiplyAt(V, index, 0, 0.0)

let rec _dot (V: array<array<float>>, i: int, accum: float) =
    // assuming vectors of same length
    if i = V.[0].Length then
        accum
    else
        let product = multiplyAt(V, i)
        _dot(V, i + 1, accum + product)

let dot (V: array<array<float>>) =
    _dot(V, 0, 0.0)

let project (u: array<float>, v: array<float>) =
    let vu = dot([|v; u|])
    let uu = dot([|u; u|])
    scale(vu / uu, u)

let rec _sumOfProjections (v: array<float>, U: list<array<float>>, i: int, accum: array<float>) =
    if i = U.Length then
        accum
    else
        let proj = project(U.[i], v)
        let sum = sum(accum, proj)
        _sumOfProjections(v, U, i + 1, sum)

let sumOfProjections (v: array<float>, U: list<array<float>>) =
    _sumOfProjections(v, U, 0, Array.zeroCreate<float> v.Length)

let magnitude(v: array<float>) =
    sqrt(dot([|v; v|]))

let normalize(v: array<float>) =
    let mag = magnitude(v)
    scale(1.0 / mag, v)

let rec _gramSchmidt (V: array<array<float>>, i: int, U: list<array<float>>) =
    if i = V.Length then
        U |> List.toArray
    elif i = 0 then
        let normalizedV = normalize(V.[i])
        _gramSchmidt(V, i + 1, U @ [normalizedV])
    else
        let pastProjections = sumOfProjections(V.[i], U)
        let diff = subtract(V.[i], pastProjections)
        let normalizedDiff = normalize(diff)
        _gramSchmidt(V, i + 1, U @ [normalizedDiff])

let validateGramSchmidtParams (V: array<array<float>>) =
    if V.Length = 0 then
        failwith "V must have at least one vector"
    else
        let length = V.[0].Length
        if V |> Array.forall (fun v -> v.Length = length) then
            ()
        else
            failwith "All vectors in V must have the same length"

let gramSchmidt (V: array<array<float>>) =
    validateGramSchmidtParams(V)
    _gramSchmidt(V, 0, [])

let main () =
    let v1 = [|1.0; 2.0; 5.0|]
    let v2 = [|3.0; 2.0; 4.0|]
    let v3 = [|3.0; 1.0; 0.0|]
    let res = gramSchmidt([|v1; v2; v3|])
    printfn "The orthonormal vectors are"
    printfn "%A" res
    let inner = dot(res)
    printfn "If the vectors are orthonormal, the inner product of the resulting vectors should be 0, or at least really close."
    printfn "The inner product of the resulting vectors is %A" inner

main()

