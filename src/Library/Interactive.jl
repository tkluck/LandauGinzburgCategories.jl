using REPL.TerminalMenus: RadioMenu, request

const WellKnownEquivalences = [
    "A₂ₙ₋₃ ≃ Dₙ (n ≥ 2) (Carqueville - Runkel 2012)" => function(n)
        TwoVariables.A{2n-3}, TwoVariables.D{n}
    end,
    "A₂ × A₂ ≃ A₅ (Recknagel - Weinreb 2017)" => function()
        TwoVariables.A₂A₂, TwoVariables.A₅
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
    "E₁₄ ≃ E₁₄ (Newton - Ros Camacho 2015)" => function()
        ThreeVariables.E₁₄{:v1}, ThreeVariables.E₁₄{:v2}
    end,
    "U₁₂ ≃ U₁₂ (Newton - Ros Camacho 2015)" => function()
        ThreeVariables.U₁₂{:v1}, ThreeVariables.U₁₂{:v3}
    end,
    "W₁₂ ≃ W₁₂ (Newton - Ros Camacho 2015)" => function()
        ThreeVariables.W₁₂{:v1}, ThreeVariables.W₁₂{:v2}
    end,
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
    return orbifold_equivalence(f, g)
end
