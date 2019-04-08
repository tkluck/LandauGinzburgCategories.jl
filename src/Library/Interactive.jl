using REPL.TerminalMenus: RadioMenu, request

const WellKnownEquivalences = [
    "A₂ₙ₋₃ ≃ Dₙ (n ≥ 2) (Carqueville - Runkel 2012)" => function(n)
        orbifold_equivalence(TwoVariables.A{2n-3}, TwoVariables.D{n})
    end,
    "A₂ × A₂ ≃ A₅ (Recknagel Weinreb 2017)" => function()
        orbifold_equivalence(TwoVariables.A₂A₂, TwoVariables.A₅)
    end,
    "A₂ × A₂ ≃ A₅ (Recknagel Weinreb 2017)" => function()
        orbifold_equivalence(TwoVariables.A₂A₂, TwoVariables.A₅)
    end,
    "D₇ ≃ E₆ (Recknagel Weinreb 2017)" => function()
        orbifold_equivalence(TwoVariables.D₇, TwoVariables.E₆)
    end,
    "D₁₀ ≃ E₇ (Recknagel Weinreb 2017)" => function()
        orbifold_equivalence(TwoVariables.D₁₀, TwoVariables.E₇)
    end,
    "D₁₆ ≃ E₈ (Recknagel Weinreb 2017)" => function()
        orbifold_equivalence(TwoVariables.D₁₆, TwoVariables.E₈)
    end,
]

function choose_equivalence()
    println("Please select one of these orbifold equivalences: ")
    options = first.(WellKnownEquivalences)
    ix = request(RadioMenu(options))

    constructor = last(WellKnownEquivalences[ix])
    method = methods(constructor).ms[1]
    if method.nargs == 1
        return constructor()
    else
        print("Please provide a value for n: ")
        n = parse(Int, readline())
        constructor(n)
    end
end
