var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#LandauGinzburgCategories.jl-Documentation-1",
    "page": "Home",
    "title": "LandauGinzburgCategories.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "This package aims to help you experiment with matrix factorizations in the bicategory of Landau-Ginzburg models. It contains implementations of operations such as unit, composition, fusion of variables, duality, and left/right quantum dimension. Future additions may include (co)evaluation maps.In addition, a library of well-known named potentials, as well as known orbifold equivalences, is provided.For getting started, have a look at Getting Started With LandauGinzburgCategories.jl."
},

{
    "location": "#Limitations-1",
    "page": "Home",
    "title": "Limitations",
    "category": "section",
    "text": "Idempotent completion. As the theory of these matrix factorizations requires an idempotent completion, not all operations can return an explicit matrix factorization. In many cases, a combination of a matrix factorization and an idempotent endomorphism is provided.Expensive computations. Many known orbifold equivalences are presented by way of an Ansatz together with the verifiable statement that the resulting equations have a solution. We represent this by placing these coefficients in a quotient ring. Every operation in this quotient ring necessarily involves a reduction step, and this can be noticably slow."
},

{
    "location": "#Table-of-contents-1",
    "page": "Home",
    "title": "Table of contents",
    "category": "section",
    "text": "Pages = [\n    # keep in sync with make.jl\n    \"index.md\",\n    \"getting-started.md\",\n    \"functions.md\",\n    \"reference.md\",\n]\nDepth = 3"
},

{
    "location": "getting-started/#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "getting-started/#Getting-Started-With-LandauGinzburgCategories.jl-1",
    "page": "Getting Started",
    "title": "Getting Started With LandauGinzburgCategories.jl",
    "category": "section",
    "text": ""
},

{
    "location": "getting-started/#Installation-1",
    "page": "Getting Started",
    "title": "Installation",
    "category": "section",
    "text": "Refer to the Julia website for details on installing Julia. As soon as you have, start it and runjulia> Pkg.add(\"https://github.com/tkluck/LandauGinzburgCategories.jl\")to install LandauGinzburgCategories and its dependencies. To test whether it worked, typeusing LandauGinzburgCategories.Library\norbifold_equivalence(TwoVariables.A{9}, TwoVariables.D{6})If you see a matrix factorization, the installation worked!"
},

{
    "location": "getting-started/#Overview-1",
    "page": "Getting Started",
    "title": "Overview",
    "category": "section",
    "text": ""
},

{
    "location": "getting-started/#Frequently-Asked-Questions-1",
    "page": "Getting Started",
    "title": "Frequently Asked Questions",
    "category": "section",
    "text": "Be the first to ask a question! Feel free to open an issue for it."
},

{
    "location": "functions/#",
    "page": "Types and Functions",
    "title": "Types and Functions",
    "category": "page",
    "text": ""
},

{
    "location": "functions/#Types-and-Functions-1",
    "page": "Types and Functions",
    "title": "Types and Functions",
    "category": "section",
    "text": ""
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.:⨷",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.:⨷",
    "category": "function",
    "text": "A⨷B\n\nGraded tensor product of ℤ/2 graded block matrices. The result consists of even/odd blocks just like the input and the operation is associative.\n\nGraded means that we introduce Koszul signs. Writing A = a^i_j e_i e^j and B = b^k_ℓ f_k f^ℓ for some (co)basis vectors, we represent the tensor product as\n\nAB = (-1)^j(k+ℓ) p^i_j q^k_ℓ e_i f_k f^ℓ e^j\n\nwhere the sign comes from commuting e^j with f_k f^ℓ. This ensures that, upon matrix multiplication between two such tensor products, the contractions between e^j and e_j, and between f^ℓ and f_ℓ, are already adjacent and do not introduce more signs.\n\nMapping the multi-indices (ik) to rows and (ℓj) to columns requires some care to ensure that we end up with even/odd blocks. A comment in the code explains how this is done. We first interlace the rows/columns to get alternating grades in rows/columns (to_alternating_grades and from_alternating_grades). Then we need only to reverse some row/column orders. Finally, we de-interlace to get even/odd blocks. (I wonder if there is a way to avoid making the detour through the alternating signs representation, but it seems hard to maintain associativity.)\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.:⨶",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.:⨶",
    "category": "function",
    "text": "A⨶B\n\nTensor product of matrix factorizations.\n\nThis is the composition operation in the bicategory or Landau Ginzburg models.\n\n\n\n\n\nQ,ϵ = ⨶(A,B,W,vars...)\n\nFinite-rank homotopy representation of AB, where we remove the variables vars from the result.\n\nWe return a matrix factorization Q together with an idempotent ϵ representing the direct summand of Q that represents AB.\n\nSee the pushforward paper by Dykerhoff&Murfet.\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.QuasiHomogeneous.centralcharge",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.QuasiHomogeneous.centralcharge",
    "category": "function",
    "text": "c = centralcharge(f, vars...)\n\nReturn the central charge of f with respect to the variables vars.\n\nThis is the value of the expression\n\nsum_i=1^n 1 - q_i\n\nwhere q_i is the grading of the ith variable under a (ℚ-valued) grading for which f is homogeneous of degree 2.\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.dual",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.dual",
    "category": "function",
    "text": "Q′ = dual(Q)\n\nReturn the dual (i.e. the adjoint) matrix factorization of Q.\n\nIf Q factors the potential f, then this is a matrix factorization that factors -f.\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.QuasiHomogeneous.find_quasihomogeneous_degrees",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.QuasiHomogeneous.find_quasihomogeneous_degrees",
    "category": "function",
    "text": "gradings = find_quasihomogeneous_degrees(f, vars...)\n\nFind the assignment of gradings to the variables vars such that f is a quasi-homogeneous polynomial with respect to these gradings, if such an assignment exists, and return it. Otherwise, raise an exception.\n\nExample\n\njulia> find_quasihomogeneous_degrees(x^4*y + x*y^9, :x, :y)\n(x = 8, y = 3)\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.fuse",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.fuse",
    "category": "function",
    "text": "docstring goes here\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.fuse_abstract",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.fuse_abstract",
    "category": "function",
    "text": "docstring goes here\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.unit_matrix_factorization",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.unit_matrix_factorization",
    "category": "function",
    "text": "unit_matrix_factorization(f; source_to_target...)\n\nA ℤ/2-graded matrix that squares to f(;source_to_target...) - f times the identity matrix.\n\nThe source for this formulation is\n\nAdjunctions and defects in Landau-Ginzburg models, Nils Carqueville and Daniel Murfet\n\n\n\n\n\n"
},

{
    "location": "functions/#Operations-1",
    "page": "Types and Functions",
    "title": "Operations",
    "category": "section",
    "text": "⨷\n⨶\ncentralcharge\ndual\nfind_quasihomogeneous_degrees\nfuse\nfuse_abstract\nunit_matrix_factorization"
},

{
    "location": "functions/#LandauGinzburgCategories.Library.orbifold_equivalence",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Library.orbifold_equivalence",
    "category": "function",
    "text": "orbifold_equivalence(f, g)\n\nReturn a matrix representing an orbifold equivalence between f and g, if one is available in the library. Return false if it is known that f and g are not orbifold equivalent. Return missing if this is not known by this library.\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Library.Potential",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Library.Potential",
    "category": "type",
    "text": "Potential{NumVars}\n\nA type representing a potential through a classification scheme. For example, the potential x^4 - y^2 is called A_3 through the ADE classification of unimodular potentials, and it is represented in Julia by the type TwoVariables.A₃. This is a subtype of Potential{2}.\n\n\n\n\n\n"
},

{
    "location": "functions/#Library-1",
    "page": "Types and Functions",
    "title": "Library",
    "category": "section",
    "text": "orbifold_equivalence\nLandauGinzburgCategories.Library.Potential"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.flatten_blocks",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.flatten_blocks",
    "category": "function",
    "text": "M = flatten_blocks(X)\n\nConstruct a matrix from a matrix of matrices by concatenating them horizontally and vertically.\n\n\n\n\n\n"
},

{
    "location": "functions/#LandauGinzburgCategories.Operations.block_diagonalization",
    "page": "Types and Functions",
    "title": "LandauGinzburgCategories.Operations.block_diagonalization",
    "category": "function",
    "text": "D, A = block_diagonalization(X)\n\nDecompose the matrix factorization X into a direct sum of irreducible matrix factorizations, and represent this direct sum as a block-diagonal matrix factorization D such that A^-1 D A = X.\n\nNOTE: for now, this function only returns D, and not yet A!\n\n\n\n\n\n"
},

{
    "location": "functions/#Internal-functions-1",
    "page": "Types and Functions",
    "title": "Internal functions",
    "category": "section",
    "text": "LandauGinzburgCategories.Operations.flatten_blocks\nLandauGinzburgCategories.Operations.block_diagonalization"
},

{
    "location": "reference/#",
    "page": "Reference Index",
    "title": "Reference Index",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Reference-Index-1",
    "page": "Reference Index",
    "title": "Reference Index",
    "category": "section",
    "text": ""
},

]}
