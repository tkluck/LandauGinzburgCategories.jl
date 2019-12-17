using REPL.TerminalMenus: RadioMenu, request

const WellKnownEquivalences = [
    "A₂ₙ₋₃ ≃ Dₙ (n ≥ 2) (Carqueville - Runkel 2012)" => function(n)
        TwoVariables.A{2n-3}, TwoVariables.D{n}
    end,
    "A₂ × A₂ ≃ A₅ (Recknagel - Weinreb 2017)" => function()
        TwoVariables.A₂A₂, TwoVariables.A₅
    end,
    "D₇ ≃ E₆ (Recknagel - Weinreb 2017)" => function()
        TwoVariables.D₇, TwoVariables.E₆
    end,
    "D₁₀ ≃ E₇ (Recknagel - Weinreb 2017)" => function()
        TwoVariables.D₁₀, TwoVariables.E₇
    end,
    # FIXME: this one is broken - something wrong with how
    # I copy-pasted the ansatz?
    #"D₁₆ ≃ E₈ (Recknagel - Weinreb 2017)" => function()
    #    TwoVariables.D₁₆, TwoVariables.E₈
    #end,
    "W₁₃ ≃ S₁₁ (Recknagel - Weinreb 2017)" => function()
        ThreeVariables.W₁₃{:v1}, ThreeVariables.S₁₁
    end,
    "E₁₃ ≃ Z₁₁ (Recknagel - Weinreb 2017)" => function()
        TwoVariables.E₁₃, TwoVariables.Z₁₁
    end,
    "Q₁₁ ≃ Z₁₃ (Recknagel - Weinreb 2017)" => function()
        ThreeVariables.Q₁₁, ThreeVariables.Z₁₃{:v1}
    end,
    "E₁₄ (v1) ≃ Q₁₀ (Newton - Ros Camacho 2015)" => function()
        ThreeVariables.E₁₄{:v1}, ThreeVariables.Q₁₀
    end,
    "E₁₄ (v2) ≃ Q₁₀ (Newton - Ros Camacho 2015)" => function()
        ThreeVariables.E₁₄{:v2}, ThreeVariables.Q₁₀
    end,
    "E₁₄ (v1) ≃ E₁₄ (v2) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.E₁₄{:v1}, ThreeVariables.E₁₄{:v2}
    end,
    "Q₁₂ (v1) ≃ Q₁₂ (v2) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.Q₁₂{:v1}, ThreeVariables.Q₁₂{:v2}
    end,
    "U₁₂ (v2) ≃ U₁₂ (v3) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.U₁₂{:v2}, ThreeVariables.U₁₂{:v3}
    end,
    "U₁₂ (v1) ≃ U₁₂ (v3) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.U₁₂{:v1}, ThreeVariables.U₁₂{:v3}
    end,
    "W₁₂ (v1) ≃ W₁₂ (v2) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.W₁₂{:v1}, ThreeVariables.W₁₂{:v2}
    end,
    "W₁₃ (v1) ≃ W₁₃ (v2) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.W₁₃{:v1}, ThreeVariables.W₁₃{:v2}
    end,
    "Z₁₃ (v1) ≃ Z₁₃ (v2) (Ros Camacho - Newton 2016)" => function()
        ThreeVariables.Z₁₃{:v1}, ThreeVariables.Z₁₃{:v2}
    end,
    "Q₁₈ ≃ E₃₀ (Kluck - Ros Camacho 2019)" => function()
        ThreeVariables.Q₁₈, ThreeVariables.E₃₀
    end,
    "E₁₈ ≃ Q₁₂ (Kluck - Ros Camacho 2019)" => function()
        ThreeVariables.E₁₈, ThreeVariables.Q₁₂
    end,
]

const Positions = [
    "Left => Right"   => (leftvars,   rightvars),
    "Left => Middle"  => (leftvars,   middlevars),
    "Middle => Right" => (middlevars, rightvars),
    "Right => Left"   => (rightvars,  leftvars),
    "Right => Middle" => (rightvars,  middlevars),
    "Middle => Left"  => (middlevars, leftvars),
]

function choose_equivalence()
    println("Please select one of these orbifold equivalences: ")
    options = first.(WellKnownEquivalences)
    ix = request(RadioMenu(options))

    closure = last(WellKnownEquivalences[ix])
    method = methods(closure).ms[1]
    if method.nargs == 1
        potentials = closure()
    else
        print("Please provide a value for n: ")
        n = parse(Int, readline())
        potentials = closure(n)
    end
    f, g = potentials

    println("Please select the how you intend to compose this equivalence: ")
    options = first.(Positions)
    ix = request(RadioMenu(options))
    leftvars, rightvars = last(Positions[ix])
    return orbifold_equivalence(f, g, leftvars(f), rightvars(g))
end

function forequivalences(fn)
    for (name, closure) in WellKnownEquivalences
        method = methods(closure).ms[1]
        if method.nargs == 1
            potentials_list = [closure()]
        else
            potentials_list = [closure(n) for n in 3:6]
        end
        for (f, g) in potentials_list
            fn(f, g)
        end
    end
end
