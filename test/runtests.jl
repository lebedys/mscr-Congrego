using DrWatson, Test
@quickactivate "Congrego"
using Congrego, JSON
using Random

include(scriptsdir("parts", "ms_original.jl"))
using .PC: pc
const external_inputs = pc.inputs

# Run test suite
println("Starting tests")
ti = time()

@testset "Congrego tests" begin
    @testset "PartCatalogue constructor" begin
        @test pc isa PartCatalogue
        @test length(pc.parts) == 35
        @test length(pc.inputs) == 2
        @test length(pc.transcription_factors) == 14

        # Minimal valid construction
        pcin = PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        @test pcin isa PartCatalogue
        @test :X ∈ pcin.transcription_factors
        @test :A ∉ pcin.transcription_factors

        # part descriptors: precomputed classification for simulate
        @test pcin.part_descriptors[:pA].class == Congrego._EXTERNALLY_ACTIVATABLE_PROMOTER
        @test pcin.part_descriptors[:pA].regulator == :A
        @test pcin.part_descriptors[:pX].class == Congrego._INTERNALLY_REPRESSIBLE_PROMOTER
        @test pcin.part_descriptors[:pX].regulator == :X
        @test pcin.part_descriptors[:X].class == Congrego._TRANSCRIPTION_FACTOR
        @test pcin.part_descriptors[:T1].class == Congrego._TERMINATOR
        @test pcin.part_descriptors[:CDS].class == Congrego._SPACER

        # verify descriptors on the full catalogue
        @test pc.part_descriptors[:pTac].class == Congrego._EXTERNALLY_ACTIVATABLE_PROMOTER
        @test pc.part_descriptors[:pTac].regulator == :IPTG
        @test pc.part_descriptors[:pTet].class == Congrego._EXTERNALLY_ACTIVATABLE_PROMOTER
        @test pc.part_descriptors[:pTet].regulator == :aTC
        @test pc.part_descriptors[:pBAD].class == Congrego._INTERNALLY_ACTIVATABLE_PROMOTER
        @test pc.part_descriptors[:pBAD].regulator == :AraC
        @test pc.part_descriptors[:pSrpR].class == Congrego._INTERNALLY_REPRESSIBLE_PROMOTER
        @test pc.part_descriptors[:pSrpR].regulator == :SrpR
        @test pc.part_descriptors[:pSrpR].is_roadblocker == true
        @test pc.part_descriptors[:pPhIF].is_roadblocker == true
        @test pc.part_descriptors[:pBM3R1].is_roadblocker == true
        @test pc.part_descriptors[:pLmrA].is_roadblocker == false
        @test pc.part_descriptors[:SrpR].class == Congrego._TRANSCRIPTION_FACTOR
        @test pc.part_descriptors[:T10].class == Congrego._TERMINATOR
        @test pc.part_descriptors[:TPhIF].class == Congrego._INDUCIBLE_TERMINATOR
        @test pc.part_descriptors[:TPhIF].regulator == :PhIF
        @test pc.part_descriptors[:tfSp].class == Congrego._SPACER

        # functional_parts: all parts except spacers and inputs
        @test pcin.functional_parts == Set([:pA, :pX, :X, :T1])
        @test :CDS ∉ pcin.functional_parts   # spacer excluded
        @test :A ∉ pcin.functional_parts      # input excluded
        @test :pTac ∈ pc.functional_parts
        @test :SrpR ∈ pc.functional_parts
        @test :T10 ∈ pc.functional_parts
        @test :TPhIF ∈ pc.functional_parts
        @test :tfSp ∉ pc.functional_parts     # spacer excluded

        # node_order: sorted (inputs ∪ TFs) then :output
        @test pcin.node_order == [:A, :X, :output]
        @test pcin.node_index == Dict(:A => 1, :X => 2, :output => 3)
        @test length(pc.node_order) == length(pc.inputs) + length(pc.transcription_factors) + 1
        @test last(pc.node_order) == :output
        @test pc.node_order[pc.node_index[:IPTG]] == :IPTG
        @test pc.node_order[pc.node_index[:SrpR]] == :SrpR
        # verify consistency: index matches position
        for (i, sym) in enumerate(pc.node_order)
            @test pc.node_index[sym] == i
        end

        ## Test assertions for valid PartCatalogue
        # 1) parts and inputs must not overlap
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X, :A]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150, :A => 50),
        )
        # 2) activatable promoters must be in parts
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:FAKE => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 3) activatable promoter transcription factors must be in parts or inputs
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :FAKE),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 4) repressible promoters must be in parts
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:FAKE => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 5) repressible promoter transcription factors must be in parts or inputs
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :FAKE),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 6) terminators must be in parts
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:FAKE]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 7) roadblockers must be in parts
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set([:FAKE]),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 8) inducible terminators must be in parts
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict(:FAKE => :X),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 9) inducible terminator transcription factors must be in parts or inputs
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict(:T1 => :FAKE),
            part_length_map = Dict(:pA => 100, :pX => 100, :CDS => 200, :T1 => 50, :X => 150),
        )
        # 10) all parts must have lengths in part_length_map
        @test_throws AssertionError PartCatalogue(
            parts = Set([:pA, :pX, :CDS, :T1, :X]),
            inputs = [:A],
            activatable_promoters = Dict(:pA => :A),
            repressible_promoters = Dict(:pX => :X),
            terminators = Set([:T1]),
            roadblockers = Set{Symbol}(),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(:pA => 100, :pX => 100),
        )
    end

    @testset "Circuit struct" begin
        @test Circuit(Symbol[]) isa Circuit
        @test length(Circuit(Symbol[])) == 0
       
        c = Circuit([:pTac, :pTet, :T10])
        @test c isa Circuit  
        @test length(c) == 3

        @test get(c, 1) == :pTac 
        @test get(c, 3) == :T10 

        io = IOBuffer()
        show(io, c)
        @test String(take!(io)) == "-- pTac -- pTet -- T10 --"
    end

    @testset "simulate returns (transcr, state)" begin
        # Single-call returns (Union{Bool,Missing}, Vector{Bool})
        # state[i] is only updated for regulatory parts (promoters, TFs);
        # spacers like tfSp fall through and stay false
        transcr, state = simulate(Circuit([:pTac, :tfSp]), :IPTG => true; pc=pc)
        @test transcr == true
        @test state isa Vector{Bool}
        @test length(state) == 2
        @test state == [true, false]  # pTac active, tfSp is a spacer (not updated)

        transcr2, state2 = simulate(Circuit([:pTac, :tfSp]), :IPTG => false; pc=pc)
        @test transcr2 == false
        @test state2 == [false, false]

        # Backward-compat Bool method also returns tuple
        transcr3, state3 = simulate(Circuit([:pTac, :tfSp]), true, false; pc=pc)
        @test transcr3 == true
        @test state3 == [true, false]

        # All-combos returns (NTuple outputs, NTuple states)
        outputs, states = simulate(Circuit([:pTac, :tfSp]), [:IPTG, :aTC]; pc=pc)
        @test outputs isa NTuple{4, Union{Bool, Missing}}
        @test states isa NTuple{4, Vector{Bool}}
        @test length(states[1]) == 2

        # String convenience also returns tuple
        outputs_s, states_s = simulate("pTac-tfSp", [:IPTG, :aTC]; pc=pc)
        @test outputs_s == outputs
        @test states_s == states

        # State is a copy (mutation-safe)
        _, s1 = simulate(Circuit([:pTac, :tfSp]), :IPTG => true; pc=pc)
        s1[1] = false
        _, s2 = simulate(Circuit([:pTac, :tfSp]), :IPTG => true; pc=pc)
        @test s2[1] == true  # fresh copy, unaffected
    end

    @testset "Circuits simulator" begin

        @test first(simulate("T10-pTac-PsrA-T10-pTac-AraC", [:IPTG, :aTC]; pc=pc)) === (false, false, true, true)
        @test first(simulate(Circuit([:T10, :pTac, :PsrA, :T10, :pTac, :AraC]), [:IPTG, :aTC]; pc=pc)) === (false, false, true, true)
        @test first(simulate("T10-pPsrA-PsrA-pTet-pBAD-tfSp", [:IPTG, :aTC]; pc=pc)) === (missing, true, missing, true)
        @test first(simulate(Circuit([:T10, :pPsrA, :PsrA, :pTet, :pBAD, :tfSp]), [:IPTG, :aTC]; pc=pc)) === (missing, true, missing, true)
        @test first(simulate("T10-pTac-PsrA-pTet-pPsrA-PsrA", [:IPTG, :aTC]; pc=pc)) === (missing, true, true, true)
        @test first(simulate(Circuit([:T10, :pTac, :PsrA, :pTet, :pPsrA, :PsrA]), [:IPTG, :aTC]; pc=pc)) === (missing, true, true, true)
        @test first(simulate("T10-pTac-PsrA-T10-pTac-tfSp", [:IPTG, :aTC]; pc=pc)) === (false, false, true, true)
        @test first(simulate(Circuit([:T10, :pTac, :PsrA, :T10, :pTac, :tfSp]), [:IPTG, :aTC]; pc=pc)) === (false, false, true, true)
        @test first(simulate("pTet-pPsrA-AraC-pTet-pPsrA-tfSp", [:IPTG, :aTC]; pc=pc)) === (true, true, true, true)
        @test first(simulate(Circuit([:pTet, :pPsrA, :AraC, :pTet, :pPsrA, :tfSp]), [:IPTG, :aTC]; pc=pc)) === (true, true, true, true)
        @test first(simulate("pTet-pPsrA-AraC-pTet-pTac-tfSp", [:IPTG, :aTC]; pc=pc)) === (true, true, true, true)
        @test first(simulate(Circuit([:pTet, :pPsrA, :AraC, :pTet, :pTac, :tfSp]), [:IPTG, :aTC]; pc=pc)) === (true, true, true, true)
        @test first(simulate("pTet-pPsrA-AraC-pTet-pPsrA-PsrA", [:IPTG, :aTC]; pc=pc)) === (missing, true, missing, true)
        @test first(simulate(Circuit([:pTet, :pPsrA, :AraC, :pTet, :pPsrA, :PsrA]), [:IPTG, :aTC]; pc=pc)) === (missing, true, missing, true)
        @test first(simulate("T10-pPsrA-PsrA-pTet-pTac-AraC", [:IPTG, :aTC]; pc=pc)) === (missing, true, true, true)
        @test first(simulate(Circuit([:T10, :pPsrA, :PsrA, :pTet, :pTac, :AraC]), [:IPTG, :aTC]; pc=pc)) === (missing, true, true, true)
        @test first(simulate("pTet-pTac-PsrA-pTet-pTac-PsrA", [:IPTG, :aTC]; pc=pc)) === (false, true, true, true)
        @test first(simulate(Circuit([:pTet, :pTac, :PsrA, :pTet, :pTac, :PsrA]), [:IPTG, :aTC]; pc=pc)) === (false, true, true, true)
        @test first(simulate("T10-pPsrA-tfSp-T10-pPsrA-tfSp", [:IPTG, :aTC]; pc=pc)) === (true, true, true, true)
        @test first(simulate(Circuit([:T10, :pPsrA, :tfSp, :T10, :pPsrA, :tfSp]), [:IPTG, :aTC]; pc=pc)) === (true, true, true, true)

        # test case with 3 cycles of oscillation (other circuits only have two cycles)
        @test first(simulate("pBM3R1-pTet-BM3R1-SrpR-pLmrA-TBM3R1-AmtR-TraRW-BM3R1", [:IPTG, :aTC]; pc=pc)) === (missing, false, missing, false)
        @test first(simulate(Circuit([:pBM3R1, :pTet, :BM3R1, :SrpR, :pLmrA, :TBM3R1, :AmtR, :TraRW, :BM3R1]), [:IPTG, :aTC]; pc=pc)) === (missing, false, missing, false)
    end

    @testset "Simulate - edge cases" begin
        # pTac requires IPTG
        @test first(simulate(Circuit([:pTac, :tfSp]), [:IPTG]; pc=pc)) == (false, true)

        # pTet requires aTC
        @test first(simulate(Circuit([:pTet, :tfSp]), [:aTC]; pc=pc)) == (false, true)

        # Terminator stops transcription
        @test first(simulate(Circuit([:T10]), false, false; pc=pc)) == false
        @test first(simulate(Circuit([:pTac, :T10, :tfSp]), true, false; pc=pc)) == false
        # New promoter starts fresh transcription after terminator
        @test first(simulate(Circuit([:pTac, :T10, :pTac, :tfSp]), true, false; pc=pc)) == true


        result, states = simulate("pTac-tfSp", [:IPTG, :aTC]; pc=pc)
        @test result isa NTuple{4, Union{Bool, Missing}}
        @test length(result) == 4
        @test states isa NTuple{4, Vector{Bool}}
    end

    @testset "Quorum sensing (pLux/AHL direct input)" begin
        # Minimal catalogue with pLux, pTac, AHL, IPTG
        pc_qs = PartCatalogue(
            parts = Set([:pTac, :pLux, :pPsrA, :PsrA, :tfSp, :T10]),
            inputs = [:IPTG, :AHL],
            activatable_promoters = Dict(:pTac => :IPTG, :pLux => :AHL),
            repressible_promoters = Dict(:pPsrA => :PsrA),
            terminators = Set([:T10]),
            roadblockers = Set([:pTac]),
            inducible_terminators = Dict{Symbol,Symbol}(),
            part_length_map = Dict(
                :pTac => 135, :pLux => 200, :pPsrA => 130,
                :PsrA => 895, :tfSp => 100, :T10 => 116,
            ),
        )

        # pLux without AHL: off
        @test first(simulate(Circuit([:pLux, :tfSp]), :AHL => false; pc=pc_qs)) == false

        # pLux with AHL: on
        @test first(simulate(Circuit([:pLux, :tfSp]), :AHL => true; pc=pc_qs)) == true

        # All combinations for pLux-tfSp (1-input: identity)
        @test first(simulate(Circuit([:pLux, :tfSp]), [:AHL]; pc=pc_qs)) == (false, true)

        # 2-input circuit: pTac drives repressor, pLux drives output
        c = Circuit([:pTac, :PsrA, :T10, :pLux, :pPsrA, :tfSp])
        # AHL=0, IPTG=0: pTac off, PsrA not made, pPsrA on, pLux off → output true (pPsrA active)
        @test first(simulate(c, :IPTG => false, :AHL => false; pc=pc_qs)) == true
        # AHL=0, IPTG=1: pTac on, PsrA made, pPsrA off, pLux off → output false
        @test first(simulate(c, :IPTG => true, :AHL => false; pc=pc_qs)) == false
        # AHL=1, IPTG=0: pTac off, PsrA not made, pPsrA on, pLux on → output true
        @test first(simulate(c, :IPTG => false, :AHL => true; pc=pc_qs)) == true
        # AHL=1, IPTG=1: pTac on, PsrA made, pPsrA off, pLux on → output true (pLux restarts)
        @test first(simulate(c, :IPTG => true, :AHL => true; pc=pc_qs)) == true

        # pLux is NOT a roadblocker — transcription continues through inactive pLux
        c_passthrough = Circuit([:pTac, :pLux, :tfSp])
        # pTac on, pLux inactive → signal passes through (not blocked)
        @test first(simulate(c_passthrough, :IPTG => true, :AHL => false; pc=pc_qs)) == true
    end

    @testset "Simulate - repressor logic" begin
        # pPsrA is ON without repressor
        @test first(simulate(Circuit([:pPsrA, :tfSp]), false, false; pc=pc)) == true

        # T10 required between repressor and repressed promoter
        # otherwise transcription continues through repressed promoter
        c = Circuit([:pTac, :PsrA, :T10, :pPsrA, :tfSp])
        @test first(simulate(c, true, false; pc=pc)) == false
        @test first(simulate(c, false, false; pc=pc)) == true
    end

    @testset "Simulate - activator logic" begin
        # pBAD requires AraC
        c = Circuit([:pTac, :AraC, :T10, :pBAD, :tfSp])
        @test first(simulate(c, false, false; pc=pc)) == false
        @test first(simulate(c, true, false; pc=pc)) == true

        # pBAD is inactive (not a roadblocker in ms_original) — signal passes when pTet is active
        c_rb = Circuit([:pTet, :pBAD, :tfSp])
        @test first(simulate(c_rb, :aTC => true; pc=pc))  == true   # pTet active, pBAD inactive and not blocking
        @test first(simulate(c_rb, :aTC => false; pc=pc)) == false  # pTet inactive, no output
    end

    @testset "Simulate - inducible terminators" begin
        # TPhIF blocks when PhIF present
        @test first(simulate(Circuit([:pTac, :PhIF, :TPhIF, :tfSp]), true, false; pc=pc)) == false
        @test first(simulate(Circuit([:pTac, :TPhIF, :tfSp]), true, false; pc=pc)) == true
    end

    @testset "output2logic" begin
        @test output2logic((false, false, false, false)) == Symbol(false)
        @test output2logic((false, false, false, true)) == :and
        @test output2logic((false, false, true, false)) == :animplyb
        @test output2logic((false, false, true, true)) == :a
        @test output2logic((false, true, false, false)) == :bnimplya
        @test output2logic((false, true, false, true)) == :b
        @test output2logic((false, true, true, false)) == :xor
        @test output2logic((false, true, true, true)) == :or
        @test output2logic((true, false, false, false)) == :nor
        @test output2logic((true, false, false, true)) == :xnor
        @test output2logic((true, false, true, false)) == :notb
        @test output2logic((true, false, true, true)) == :bimplya
        @test output2logic((true, true, false, false)) == :nota
        @test output2logic((true, true, false, true)) == :aimplyb
        @test output2logic((true, true, true, false)) == :nand
        @test output2logic((true, true, true, true)) == Symbol(true)
    end

    @testset "sequencing_price" begin
        @test sequencing_price(1) == 37.0
        @test sequencing_price(100) == 37.0
        @test sequencing_price(501) == 54.0
        @test sequencing_price(2000) == 150.0
        @test sequencing_price(3000) == 245.0
        @test sequencing_price(3001) == Inf
    end

    @testset "Simulation - stress tests" begin
        @test first(simulate(Circuit(rand(collect(pc.parts), 15)), [:IPTG, :aTC]; pc=pc)) isa NTuple{4, Union{Bool, Missing}}
        @test first(simulate(Circuit([:pTac, :AraC, :pBAD, :PsrA, :pPsrA, :T10, :pTet, :PhIF, :TPhIF, :tfSp]), [:IPTG, :aTC]; pc=pc)) isa NTuple{4, Union{Bool, Missing}}
        @test first(simulate(Circuit([:pTac, :PsrA, :pPsrA, :PhIF, :pPhIF, :BM3R1, :pBM3R1, :tfSp]), [:IPTG, :aTC]; pc=pc)) isa NTuple{4, Union{Bool, Missing}}
    end

    @testset "Simulate - v1 backward compatibility" begin
        # Verify refactored simulate matches archived _simulate_v1 for all test circuits
        test_circuits = [
            Circuit([:pTac, :tfSp]),
            Circuit([:pTet, :tfSp]),
            Circuit([:T10]),
            Circuit([:pTac, :T10, :tfSp]),
            Circuit([:pTac, :T10, :pTac, :tfSp]),
            Circuit([:pPsrA, :tfSp]),
            Circuit([:pTac, :PsrA, :T10, :pPsrA, :tfSp]),
            Circuit([:pTac, :AraC, :T10, :pBAD, :tfSp]),
            Circuit([:pTac, :PhIF, :TPhIF, :tfSp]),
            Circuit([:pTac, :TPhIF, :tfSp]),
            Circuit(Symbol.(split("T10-pTac-PsrA-T10-pTac-AraC", "-"))),
            Circuit(Symbol.(split("pTet-pPsrA-AraC-pTet-pPsrA-PsrA", "-"))),
            Circuit(Symbol.(split("pTet-pTac-PsrA-pTet-pTac-PsrA", "-"))),
            Circuit(Symbol.(split("T10-pPsrA-tfSp-T10-pPsrA-tfSp", "-"))),
        ]
        for c in test_circuits
            for IPTG in [false, true], aTC in [false, true]
                @test first(simulate(c, IPTG, aTC; pc=pc)) === Congrego._simulate_v1(c, IPTG, aTC; pc=pc)
            end
        end
    end

    @testset "Simulate - Pair interface" begin
        c = Circuit([:pTac, :tfSp])
        @test simulate(c, :IPTG => true, :aTC => false; pc=pc) == simulate(c, true, false; pc=pc)
        @test simulate(c, :IPTG => false, :aTC => false; pc=pc) == simulate(c, false, false; pc=pc)
        @test simulate(c, :IPTG => false, :aTC => true; pc=pc) == simulate(c, false, true; pc=pc)
        # Omitted inputs default to false
        @test simulate(c, :IPTG => true; pc=pc) == simulate(c, true, false; pc=pc)
        @test simulate(c, :aTC => true; pc=pc) == simulate(c, false, true; pc=pc)
        # Order-independent
        @test simulate(c, :aTC => false, :IPTG => true; pc=pc) == simulate(c, :IPTG => true, :aTC => false; pc=pc)
        # All-combos with explicit input list
        @test first(simulate(c, [:IPTG, :aTC]; pc=pc)) isa NTuple{4, Union{Bool, Missing}}
    end

    @testset "external_inputs" begin
        @test :IPTG ∈ external_inputs
        @test :aTC ∈ external_inputs
        # No external inputs should appear in the parts set:
        @test isempty(Set(external_inputs) ∩ pc.parts)
    end

    @testset "n_phenotypes" begin
        @test n_phenotypes(1) == 4
        @test n_phenotypes(2) == 16
        @test n_phenotypes(3) == 256
    end

    @testset "pheno_index" begin
        @test pheno_index((false, false, false, false)) == 1
        @test pheno_index((false, false, false, true)) == 2
        @test pheno_index((true, true, false, false)) == 13
        # Must match output2logic indexing for all 2-input cases to maintain backwards-compatibility
        for combo in Iterators.product(fill([false,true], 4)...)
            @test pheno_index(Tuple(combo)) == findfirst(==(output2logic(Tuple(combo))), logic_functions)
        end
    end

    @testset "pheno_label" begin
        # 1-input
        @test pheno_label((false, false)) == Symbol(false)
        @test pheno_label((false, true)) == :identity
        @test pheno_label((true, false)) == :not
        @test pheno_label((true, true)) == Symbol(true)
        # 2-input: delegates to output2logic
        @test pheno_label((false, false, false, true)) == :and
        @test pheno_label((true, true, true, false)) == :nand
        # 3-input: returns Symbol of binary string
        @test pheno_label((false, false, false, false, true, true, true, true)) == Symbol("00001111")
        @test pheno_label((true, true, true, true, true, true, true, true)) == Symbol("11111111")
    end

    @testset "characterize" begin
        # Produces 2 simple circuits: pTac-tfSp and pTet-tfSp
        g = Generator([
            [:pTac, :pTet],
            [:tfSp]
        ])
        counts, n = characterize(g, [:IPTG, :aTC]; pc=pc)
        @test n == 2
        @test length(counts) == n_phenotypes(2)
        @test sum(counts) == n  # no oscillators
        @test counts[pheno_index((false, false, true, true))] == 1  # :a (pTac-tfSp)
        @test counts[pheno_index((false, true, false, true))] == 1  # :b (pTet-tfSp)


        # Produces 4 circuits, 2 of which oscillate
        g = Generator([
            [:pTet, :pTac],
            [:pPsrA, :tfSp],
            [:PsrA]
        ])
        counts, n = characterize(g, [:IPTG, :aTC]; pc=pc)
        @test n == 4
        @test length(counts) == n_phenotypes(2)
        # oscillators don't get counted
        @test sum(counts) == n - 2
        @test counts[pheno_index((false, false, true, true))] == 1  # :a (pTac-tfSp)
        @test counts[pheno_index((false, true, false, true))] == 1  # :b (pTet-tfSp)


        # Reference parts used to build Michiel's library.
        # Taken from 'final_ordered_generator.txt'
        g = Generator([
            [[:pTac, :SrpR]],
            [[:pTet, :pTac, :TSrpR, :BM3R1],
             [:pTet, :pTac, :TSrpR],
             [:pTet]],
            Vector{Symbol}[],
            Vector{Symbol}[],
            Vector{Symbol}[],
            [[:LmrA], [:pTra], [:tfSp], [:TraRW]],
            [[:AraC, :PhIF]],
            [[:pSrpR, :pTra],
             [:pTac, :TSrpR],
             [:BetI],
             [:pSrpR, :pBAD]],
            [[:BetI, :pTet, :BM3R1, :IcaRA, :pBM3R1, :pBAD],
             [:BetI, :pTet, :TPhIF],
             [:BetI, :pTet, :BM3R1, :IcaRA, :pBM3R1],
             [:BM3R1],
             [:pTac, :TPhIF],
             [:BetI, :pTet]],
            Vector{Symbol}[],
        ])
        # Capable of generating 288 circuits
        # Building cost: 2880.0 GBP
        @test number_of_circuits(g) == 288
        @test get_price(g, pc) == 1552.0

        #             ┌                                        ┐ 
        #       false ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 21   
        #         and ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #    animplyb ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 19      
        #           a ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #    bnimplya ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 20     
        #           b ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 20     
        #         xor ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #          or ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #         nor ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 20     
        #        xnor ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #        notb ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 20     
        #    animplyb ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #        nota ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 19      
        #    bnimplya ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 19      
        #        nand ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■ 16            
        #        true ┤■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 18        
        #             └                                        ┘ 
        counts, n = characterize(g, [:IPTG, :aTC]; pc=pc)
        @test n == 288
        @test sum(counts) == 288
        @test counts == [21, 16, 19, 16, 20, 20, 16, 16, 20, 16, 20, 16, 19, 19, 16, 18]
    end

    @testset "Generator utilities" begin
        # Empty generator
        g_empty = Generator([Symbol[], Symbol[]])
        @test isempty(g_empty)
        @test number_of_slots(g_empty) == 0
        @test number_of_parts(g_empty) == 0

        # CGenerator with all-empty slots — iter must return Vector{Symbol}
        g_empty_composite = Generator([Vector{Symbol}[], Vector{Symbol}[]])
        @test Congrego.part_type(g_empty_composite) == Congrego.CompositeParts
        for c in Congrego.iter(g_empty_composite)
            @test c isa Vector{Symbol}
        end

        # Simple generator with no composite parts
        g = Generator([[:pTac, :pTet], [:tfSp], Symbol[]])

        @test number_of_slots(g) == 2  # empty slot ignored
        @test number_of_parts(g) == 3  # empty slot ignored
        @test number_of_circuits(g) == 2
        @test !isempty(g)


        # gen2parts: collect parts from each slot
        gen_parts = gen2parts(g)
        @test length(gen_parts) == 3  # empty slot included
        @test Set(gen_parts[1]) == Set([:pTac, :pTet])
        @test Set(gen_parts[2]) == Set([:tfSp])
        @test isempty(gen_parts[3])


        # get_parts: unique sorted parts across all slots
        unique_gen_parts = get_parts(g)
        @test issorted(unique_gen_parts)
        @test length(unique_gen_parts) == length(unique(unique_gen_parts))
        @test Set(unique_gen_parts) == Set([:pTac, :pTet, :tfSp])



        # Generator with composite parts
        g_composite = Generator([
            [[:pTac, :SrpR]],  # pTac-SrpR composite
            [[:pTet]],
            Vector{Symbol}[],
        ])

        @test number_of_slots(g_composite) == 2
        @test number_of_parts(g_composite) == 2

        # gen2parts: collect composite parts from each slot
        composite_gen_parts = gen2parts(g_composite)
        @test length(composite_gen_parts) == 3  # includes empty slot
        @test [:pTac, :SrpR] ∈ composite_gen_parts[1]
        @test [:pTet] ∈ composite_gen_parts[2]
        @test isempty(composite_gen_parts[3])

        # get_parts for composite: flattens all sub-parts
        unique_composite_gen_parts = get_parts(g_composite)
        @test issorted(unique_composite_gen_parts)
        @test Set(unique_composite_gen_parts) == Set([:pTac, :pTet, :SrpR])


        # price comparison
        g_small = Generator([[:pTac], [:tfSp]])
        g_big = Generator([[:pTac, :pTet], [:tfSp, :AraC]])
        @test get_price(g_big, pc) > get_price(g_small, pc)
    end

    @testset "Algorithms" begin
        @testset "tournament_selection" begin
            pop = [Generator([[:pTac], [:tfSp]]),
                   Generator([[:pTet], [:tfSp]]),
                   Generator([[:pBAD], [:tfSp]]),
                   Generator([[:pTra], [:tfSp]])]
            objectives = [1.0, 3.0, 2.0, 0.5]

            # check return size
            @test length(tournament_selection(objectives, pop)) == length(pop)
            @test length(tournament_selection(objectives, pop; elitism=false)) == length(pop)
        end

        @testset "crossover" begin
            g1 = Generator([[:pTac, :pTet], [:tfSp]])
            g2 = Generator([[:pBAD], [:AraC]])
            c1, c2 = crossover(g1, g2)

            # slot number preserved
            @test length(c1.slots) == length(g1.slots)
            @test length(c2.slots) == length(g2.slots)

            # original generators are not mutated (deepcopy inside)
            @test g1.slots[1] == Set([:pTac, :pTet])
            @test g2.slots[1] == Set([:pBAD])
        end

        @testset "mutate - no composite parts" begin
            Random.seed!(42)
            parts = collect(pc.parts)
            g = Generator([[:pTac, :pTet], [:tfSp], Symbol[]])
            @test Congrego.part_type(g) == Congrego.SingletonParts

            # with pmut=0.0 nothing changes
            g_nomut = mutate(g, parts; pmut=0.0)
            @test g_nomut.slots == g.slots

            # with pmut=1.0 every slot is mutated
            g_mut = mutate(g, parts; pmut=1.0)
            @test g_mut isa Generator
            @test length(g_mut.slots) == length(g.slots)

            # original not mutated
            @test g.slots == Generator([[:pTac, :pTet], [:tfSp], Symbol[]]).slots
        end

        @testset "mutate - composite parts" begin
            Random.seed!(42)
            parts = collect(pc.parts)
            g = Generator([
                [[:pTac, :SrpR]],
                [[:pTet]],
                Vector{Symbol}[]
            ])
            @test Congrego.part_type(g) == Congrego.CompositeParts

            # with pmut=0.0 nothing changes
            g_nomut = mutate(g, parts; pmut=0.0)
            @test g_nomut.slots == g.slots

            # with pmut=1.0 every slot is mutated
            g_mut = mutate(g, parts; pmut=1.0)
            @test g_mut isa Generator
            @test length(g_mut.slots) == length(g.slots)

            # original not mutated
            @test g.slots == Generator([
                [[:pTac, :SrpR]],
                [[:pTet]],
                Vector{Symbol}[]
            ]).slots
        end

        @testset "mutate! - singleton parts (in-place)" begin
            Random.seed!(42)
            parts = collect(pc.parts)
            g = Generator([[:pTac, :pTet], [:tfSp], Symbol[]])

            # mutate! returns the same object (not a copy)
            g_mut = mutate!(g, parts; pmut=1.0)
            @test g_mut === g
        end

        @testset "mutate! - composite parts (in-place)" begin
            Random.seed!(42)
            parts = collect(pc.parts)
            g = Generator([
                [[:pTac, :SrpR]],
                [[:pTet]],
                Vector{Symbol}[]
            ])

            g_mut = mutate!(g, parts; pmut=1.0)
            @test g_mut === g
        end

        @testset "recombinate!" begin
            parts = collect(pc.parts)
            pop = [Generator([[:pTac], [:tfSp]]),
                   Generator([[:pTet], [:tfSp]]),
                   Generator([[:pBAD], [:AraC]]),
                   Generator([[:pTra], [:tfSp]])]

            result = recombinate!(copy(pop), parts; pmut=0.05, pcrossover=0.25)

            # population size preserved
            @test length(result) == length(pop)
            # all elements are Generators
            @test all(g -> g isa Generator, result)
        end

        @testset "genetic_algorithm! smoke" begin
            Random.seed!(123)
            parts = collect(pc.parts)
            pop = [Generator([[:pTac, :pTet], [:tfSp]]) for _ in 1:8]

            fitness(g) = Float64(number_of_circuits(g))

            pop_out, objectives, best_obj = genetic_algorithm!(
                fitness, pop, parts;
                n_gen=5, pmut=0.1, pcrossover=0.25, verbose=false)

            @test length(pop_out) == 8
            @test all(g -> g isa Generator, pop_out)
            @test length(objectives) == 8
            @test all(x -> x ≥ 0, objectives)
            @test length(best_obj) == 6  # n_gen + 1
            @test all(x -> x ≥ 0, best_obj)
        end

        @testset "map_elites smoke" begin
            Random.seed!(456)
            parts = collect(pc.parts)
            pop = [Generator([[:pTac, :pTet], [:tfSp]]) for _ in 1:8]

            fitness(g) = (Float64(number_of_circuits(g)), (0.0, 0.0, 0.0, 0.0))
            descriptors(g) = (number_of_slots(g), number_of_parts(g))

            archive, best_obj, history = map_elites(
                fitness, pop, descriptors, parts;
                n_iter=5, pmut=0.1, pcrossover=0.25, verbose=false)

            @test archive isa Dict
            @test !isempty(archive)
            # archive values are (fitness, Generator, components) tuples
            for (k, (f, g, c)) in archive
                @test f isa Float64
                @test f ≥ 0
                @test g isa Generator
                @test c isa NTuple{4, Float64}
                @test k isa Tuple{Int, Int}
            end
            @test length(best_obj) == 6  # n_iter + 1
            @test all(x -> x ≥ 0, best_obj)
            # history: each cell has at least one entry
            @test history isa Dict
            @test Set(keys(history)) == Set(keys(archive))
            for (k, entries) in history
                @test length(entries) ≥ 1
                # entries are sorted by iteration and fitness is non-decreasing
                for (iter, fv, comp) in entries
                    @test iter isa Int
                    @test fv isa Float64
                    @test comp isa NTuple{4, Float64}
                end
            end
        end
    end

    @testset "PartCatalogue JSON serialisation" begin
        d = JSON.lower(pc)

        # all expected keys present
        @test Set(string.(keys(d))) == Set([
            "parts", 
            "inputs", 
            "activatable_promoters", "repressible_promoters",
            "terminators", "roadblockers", "inducible_terminators",
            "promoters", "transcription_factors", 
            "part_length_map"
        ])

        # determinism — serialise twice, compare
        @test JSON.lower(pc) == JSON.lower(pc)

        # pc → JSON → PartCatalogue(::Dict) equivalence check
        json_str = JSON.json(d)
        d2 = JSON.parse(json_str)
        pc2 = PartCatalogue(d2)
        @test JSON.lower(pc2) == d
    end

    @testset "Levenshtein edit distance" begin
        verbose=true

        # empty
        @test levenshtein(Symbol[], Symbol[]; verbose) == 0
        @test normalised_levenshtein(Symbol[], Symbol[]; verbose) == 0.0
        
        # no difference
        @test levenshtein([:a, :b, :c], [:a, :b, :c]; verbose) == 0
        @test normalised_levenshtein([:a, :b, :c], [:a, :b, :c]; verbose) == 0.0

        # full substitution
        @test levenshtein([:a, :b, :c], [:x, :y, :z]; verbose) == 3
        @test normalised_levenshtein([:a, :b, :c], [:x, :y, :z]; verbose) == 1.0
        
        # single addition
        @test levenshtein([:a, :b], [:a, :b, :c]; verbose) == 1
        @test normalised_levenshtein([:a, :b], [:a, :b, :c]; verbose) ≈ 1/3
        
        # single substitution
        @test levenshtein([:a, :b, :c], [:a, :y, :c]; verbose) == 1
        @test normalised_levenshtein([:a, :b, :c], [:a, :y, :c]; verbose) ≈ 1/3
        
        # single addition vs empty
        @test levenshtein([:a], Symbol[]; verbose) == 1
        @test normalised_levenshtein([:a], Symbol[]; verbose) == 1.0

        # verbose mode doesn't error
        @test levenshtein([:a, :b], [:a, :c]; verbose) == 1

        verbose=false
    end

    @testset "characterize_with_diversity" begin
        "Test edge cases for archiectural diversity"
        verbose=false

        # No diversity in a single generator
        g_single = Generator([[:pTac], [:tfSp]])
        counts_d, n_d, arch_div_d = characterize_with_diversity(g_single, [:IPTG, :aTC]; pc=pc)
        counts_ref, n_ref = characterize(g_single, [:IPTG, :aTC]; pc=pc)
        @test counts_d == counts_ref
        @test n_d == n_ref
        @test arch_div_d == 0.0

        # No diversity between two circuits in different phenotypes
        g_diverse = Generator([[:pTac, :pTet], [:tfSp]])
        counts_d2, n_d2, arch_div_d2 = characterize_with_diversity(g_diverse, [:IPTG, :aTC]; pc=pc)
        counts_ref2, n_ref2 = characterize(g_diverse, [:IPTG, :aTC]; pc=pc)
        @test counts_d2 == counts_ref2
        @test n_d2 == n_ref2
        @test n_d2 == 2
        @test arch_div_d2 == 0.0

        g_mixed = Generator([
            [:pTet, :pTac],
            [:pPsrA, :tfSp],
            [:PsrA]
        ])
        counts_d3, n_d3, arch_div_d3 = characterize_with_diversity(g_mixed, [:IPTG, :aTC]; pc=pc)
        counts_ref3, n_ref3 = characterize(g_mixed, [:IPTG, :aTC]; pc=pc)
        @test counts_d3 == counts_ref3
        @test n_d3 == n_ref3
        @test n_d3 == 4

        # Works for composite parts
        g_comp = Generator([
            [[:pTac, :SrpR]],
            [[:pTet]],
            Vector{Symbol}[],
        ])
        counts_d4, n_d4, arch_div_d4 = characterize_with_diversity(g_comp, [:IPTG, :aTC]; pc=pc)
        counts_ref4, n_ref4 = characterize(g_comp, [:IPTG, :aTC]; pc=pc)
        @test counts_d4 == counts_ref4
        @test n_d4 == n_ref4
        # TODO - log the arch_div value when verbose
        @test 0.0 ≤ arch_div_d4 ≤ 1.0

        # TODO - make this bigger
        # arch_div ϵ [0, 1]
        g_big = Generator([
            [:pTac, :pTet],
            [:tfSp, :PsrA],
            [:T10, :pPsrA]
        ])
        _, _, arch_div_big = characterize_with_diversity(g_big, [:IPTG, :aTC]; pc=pc)
        @test 0.0 ≤ arch_div_big ≤ 1.0

        # ── metric kwarg ─────────────────────────────────────────────
        # default (:edit_distance) matches existing behavior
        counts_ed, n_ed, div_ed = characterize_with_diversity(g_big, [:IPTG, :aTC]; pc=pc, metric=:edit_distance)
        @test div_ed == arch_div_big    # explicit :edit_distance == default

        # :graph_distance returns valid diversity in [0, 1]
        counts_gd, n_gd, div_gd = characterize_with_diversity(g_big, [:IPTG, :aTC]; pc=pc, metric=:graph_distance)
        @test counts_gd == counts_ed    # phenotype counts unchanged by metric
        @test n_gd == n_ed
        @test 0.0 ≤ div_gd ≤ 1.0

        # single-circuit generator → diversity 0 for both metrics
        _, _, div_single_ed = characterize_with_diversity(g_single, [:IPTG, :aTC]; pc=pc, metric=:edit_distance)
        _, _, div_single_gd = characterize_with_diversity(g_single, [:IPTG, :aTC]; pc=pc, metric=:graph_distance)
        @test div_single_ed == 0.0
        @test div_single_gd == 0.0

        # Subsampling path: force n > max_sample. The constant-false phenotype
        # bin here holds 7 of the 8 circuits (every arrangement with ≥1 T10),
        # so max_sample=2 exercises _get_subsampled_indices + the pairwise loop.
        g_sub = Generator([[:pTac], [:tfSp, :T10], [:tfSp, :T10], [:tfSp, :T10]])
        for metric in (:edit_distance, :graph_distance)
            _, _, d_sub1 = characterize_with_diversity(g_sub, [:IPTG, :aTC]; pc=pc, max_sample=2, metric=metric)
            _, _, d_sub2 = characterize_with_diversity(g_sub, [:IPTG, :aTC]; pc=pc, max_sample=2, metric=metric)
            @test 0.0 ≤ d_sub1 ≤ 1.0
            @test d_sub1 == d_sub2   # deterministic
        end
        verbose=false
    end

    @testset "_get_subsampled_indices" begin
        # n ≤ max_sample → identity passthrough
        @test Congrego._get_subsampled_indices(10, 50) == collect(1:10)
        @test Congrego._get_subsampled_indices(50, 50) == collect(1:50)

        # n > max_sample → exactly max_sample evenly-spaced indices spanning the
        # full range incl. both endpoints. Regression for the old cld-stride bug
        # that under-sampled (76 not 100) and could duplicate index n.
        idx = Congrego._get_subsampled_indices(301, 100)
        @test length(idx) == 100
        @test allunique(idx)
        @test issorted(idx)            # issorted + allunique ⟹ strictly increasing
        @test first(idx) == 1
        @test last(idx) == 301

        # edge cases (hit the existing length<1 and length==1 guards upstream)
        @test Congrego._get_subsampled_indices(0, 50) == Int[]
        @test Congrego._get_subsampled_indices(1, 50) == [1]
    end

    @testset "normalised_graph_distance" begin
        # hand-crafted 3×3 matrices for focused distance testing
        A = zeros(Int8, 3, 3)
        B = zeros(Int8, 3, 3)

        # identical zero matrices → 0.0
        @test Congrego.normalised_graph_distance(A, B) == 0.0

        # identical non-zero matrices → 0.0
        A[1,2] = Int8(1); A[2,3] = Int8(-1)
        B[1,2] = Int8(1); B[2,3] = Int8(-1)
        @test Congrego.normalised_graph_distance(A, B) == 0.0

        # one edge differs in sign → 1 mismatch / 2 active = 0.5
        B[2,3] = Int8(1)  # was -1, now +1
        @test Congrego.normalised_graph_distance(A, B) ≈ 0.5

        # completely disjoint edges → 1.0
        fill!(A, Int8(0)); fill!(B, Int8(0))
        A[1,2] = Int8(1)
        B[2,3] = Int8(-1)
        @test Congrego.normalised_graph_distance(A, B) == 1.0

        # edge present in A only → 1 mismatch / 1 active = 1.0
        fill!(A, Int8(0)); fill!(B, Int8(0))
        A[1,2] = Int8(1)
        @test Congrego.normalised_graph_distance(A, B) == 1.0

        # symmetry
        fill!(A, Int8(0)); fill!(B, Int8(0))
        A[1,2] = Int8(1); A[1,3] = Int8(-1)
        B[1,2] = Int8(1); B[2,3] = Int8(1)
        @test Congrego.normalised_graph_distance(A, B) == Congrego.normalised_graph_distance(B, A)
    end

    @testset "part_length_bp" begin
        @test part_length_bp(:pTac, pc) == pc.part_length_map[:pTac]
        @test part_length_bp([:pTac, :SrpR], pc) == pc.part_length_map[:pTac] + pc.part_length_map[:SrpR]
    end

    @testset "mean_part_bp" begin
        # Singleton: mean of individual part lengths
        g1 = Generator([[:pTac, :pTet], [:tfSp]])
        expected = (pc.part_length_map[:pTac] + pc.part_length_map[:pTet] + pc.part_length_map[:tfSp]) / 3
        @test mean_part_bp(g1, pc) ≈ expected

        # Composite: mean of composite part lengths
        g2 = Generator([[[:pTac, :SrpR], [:pTet]], [[:tfSp]]])
        len_a = pc.part_length_map[:pTac] + pc.part_length_map[:SrpR]
        len_b = pc.part_length_map[:pTet]
        len_c = pc.part_length_map[:tfSp]
        @test mean_part_bp(g2, pc) ≈ (len_a + len_b + len_c) / 3

        # Empty generator
        g3 = Generator([Symbol[], Symbol[]])
        @test mean_part_bp(g3, pc) == 0.0
    end

    @testset "get_price refactored" begin
        @test get_price(:pTac, pc) == sequencing_price(pc.part_length_map[:pTac])
        @test get_price([:pTac, :SrpR], pc) == sequencing_price(pc.part_length_map[:pTac] + pc.part_length_map[:SrpR])
    end

    @testset "utils" begin
        @testset "best_cheapest_generator" begin
            g_cheap     = Generator([[:tfSp]])            # :tfSp = 100 bp             → £67
            g_medium    = Generator([[:IcaRA]])           # :IcaRA = 743 bp            → £101
            g_expensive = Generator([[:BM3R1]])           # :BM3R1 = 751 bp            → £160
            g_composite = Generator([[[:IcaRA, :BetI]]])  # composite: 743+749=1492 bp → £206
            @assert get_price(g_cheap,     pc) < get_price(g_medium,    pc)
            @assert get_price(g_medium,    pc) < get_price(g_expensive, pc)
            @assert get_price(g_expensive, pc) < get_price(g_composite, pc)

            dummy_c = (0.0, 0.0, 0.0, 0.0)

            # Case 1: winner by fitness (cost is irrelevant)
            archive = Dict(
                (1, 250)  => (15.0, g_composite, dummy_c), # fitness wins
                (2, 500)  => (10.0, g_expensive, dummy_c),
                (3, 750)  => ( 5.0, g_medium, dummy_c),
                (4, 1000) => ( 1.0, g_cheap, dummy_c),
            )
            f, g, _ = best_cheapest_generator(archive, pc)
            @test f == 15.0
            @test g === g_composite

            # Case 2: tied fitness, cheapest wins
            archive_tied = Dict(
                (1, 250)  => (10.0, g_expensive, dummy_c),
                (2, 500)  => (10.0, g_medium, dummy_c),
                (3, 750)  => (10.0, g_cheap, dummy_c), # lowest cost wins
                (4, 1000) => ( 5.0, g_composite, dummy_c),
            )
            f2, g2, _ = best_cheapest_generator(archive_tied, pc)
            @test f2 == 10.0
            @test g2 === g_cheap
        end
    end

    # ── extract_adjacency! ─────────────────────────────────────────────
    # Load qs2 catalogue (pTac is a roadblocker here, unlike ms_original)
    include(scriptsdir("parts", "qs2.jl"))
    pc_qs2 = PC.pc

    @testset "extract_adjacency!" begin
        N = length(pc_qs2.node_order)
        adj = zeros(Int8, N, N)
        ni = pc_qs2.node_index

        # ── Group 1: Basic single-class ──────────────────────────────
        @testset "spacer only" begin
            Congrego.extract_adjacency!(adj, [:tfSp]; pc=pc_qs2)
            @test all(adj .== 0)
        end

        @testset "single promoter → output edge" begin
            Congrego.extract_adjacency!(adj, [:pTac]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test count(!=(0), adj) == 1
        end

        @testset "promoter + TF" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test count(!=(0), adj) == 2
        end

        @testset "repressible promoter + TF" begin
            Congrego.extract_adjacency!(adj, [:pSrpR, :BM3R1]; pc=pc_qs2)
            @test adj[ni[:SrpR], ni[:BM3R1]] == Int8(-1)
            @test adj[ni[:SrpR], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 2
        end

        @testset "TF with no promoter in scope" begin
            Congrego.extract_adjacency!(adj, [:SrpR]; pc=pc_qs2)
            @test all(adj .== 0)
        end

        @testset "terminator with nothing to terminate" begin
            Congrego.extract_adjacency!(adj, [:T10]; pc=pc_qs2)
            @test all(adj .== 0)
        end

        # ── Group 2: Roadblocker vs non-roadblocker ──────────────────
        @testset "RB then non-RB: both contribute (OR)" begin
            # pTac=RB, pBetI=non-RB → both in scope
            Congrego.extract_adjacency!(adj, [:pTac, :pBetI, :SrpR]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:BetI], ni[:SrpR]] == Int8(-1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test adj[ni[:BetI], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 4
        end

        @testset "non-RB then RB: roadblocker replaces" begin
            # pBetI=non-RB, pTac=RB → pTac replaces pBetI's scope
            Congrego.extract_adjacency!(adj, [:pBetI, :pTac, :SrpR]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test count(!=(0), adj) == 2
        end

        @testset "two non-RBs: both contribute (OR)" begin
            # pBetI=non-RB, pLmrA=non-RB → both in scope
            Congrego.extract_adjacency!(adj, [:pBetI, :pLmrA, :SrpR]; pc=pc_qs2)
            @test adj[ni[:BetI], ni[:SrpR]] == Int8(-1)
            @test adj[ni[:LmrA], ni[:SrpR]] == Int8(-1)
            @test adj[ni[:BetI], ni[:output]] == Int8(-1)
            @test adj[ni[:LmrA], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 4
        end

        # ── Group 3: Terminators ─────────────────────────────────────
        @testset "terminator separates two segments" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :T10, :pSrpR, :BM3R1]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:SrpR], ni[:BM3R1]] == Int8(-1)
            @test adj[ni[:SrpR], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 3
        end

        @testset "terminator kills scope → TF unreachable" begin
            Congrego.extract_adjacency!(adj, [:pTac, :T10, :SrpR]; pc=pc_qs2)
            @test all(adj .== 0)
        end

        @testset "terminator at end → no output edges" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :T10]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test count(!=(0), adj) == 1
        end

        # ── Group 4: Inducible terminators ───────────────────────────
        @testset "inducible terminator adds repressive edge" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :TBM3R1, :BetI]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:IPTG], ni[:BetI]] == Int8(1)
            @test adj[ni[:BM3R1], ni[:BetI]] == Int8(-1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test adj[ni[:BM3R1], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 5
        end

        @testset "contrast: T10 clears vs inducible preserves" begin
            # T10 clears → BetI unreachable
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :T10, :BetI]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test count(!=(0), adj) == 1
        end

        @testset "inducible terminator before any promoter → guard" begin
            # TSrpR at start: scope empty → guard prevents spurious edge
            Congrego.extract_adjacency!(adj, [:TSrpR, :pTac, :SrpR]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test count(!=(0), adj) == 2
        end

        @testset "T10 then inducible terminator → guard fires" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :T10, :TSrpR, :BM3R1]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test count(!=(0), adj) == 1
        end

        @testset "two consecutive inducible terminators accumulate" begin
            Congrego.extract_adjacency!(adj, [:pTac, :TSrpR, :TBM3R1, :BetI]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:BetI]] == Int8(1)
            @test adj[ni[:SrpR], ni[:BetI]] == Int8(-1)
            @test adj[ni[:BM3R1], ni[:BetI]] == Int8(-1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test adj[ni[:SrpR], ni[:output]] == Int8(-1)
            @test adj[ni[:BM3R1], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 6
        end

        @testset "RB promoter clears inducible terminator contribution" begin
            # pBM3R1 is RB → empties scope including TSrpR's entry
            Congrego.extract_adjacency!(adj, [:pTac, :TSrpR, :pBM3R1, :BetI]; pc=pc_qs2)
            @test adj[ni[:BM3R1], ni[:BetI]] == Int8(-1)
            @test adj[ni[:BM3R1], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 2
        end

        # ── Group 5: Shared regulators ───────────────────────────────
        @testset "pSrpR + TSrpR: duplicate write is idempotent" begin
            Congrego.extract_adjacency!(adj, [:pSrpR, :TSrpR, :BM3R1]; pc=pc_qs2)
            @test adj[ni[:SrpR], ni[:BM3R1]] == Int8(-1)
            @test adj[ni[:SrpR], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 2
        end

        @testset "cascade then shared regulator" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :pSrpR, :TSrpR, :BM3R1]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:SrpR], ni[:BM3R1]] == Int8(-1)
            @test adj[ni[:SrpR], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 3
        end

        # ── Group 6: Complex multi-segment circuits ──────────────────
        @testset "standard cascade" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :pSrpR, :BM3R1, :pBM3R1, :tfSp]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:SrpR], ni[:BM3R1]] == Int8(-1)
            @test adj[ni[:BM3R1], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 3
        end

        @testset "reordered cascade → different edges" begin
            Congrego.extract_adjacency!(adj, [:pTac, :BM3R1, :pBM3R1, :SrpR, :pSrpR, :tfSp]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:BM3R1]] == Int8(1)
            @test adj[ni[:BM3R1], ni[:SrpR]] == Int8(-1)
            @test adj[ni[:SrpR], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 3
        end

        @testset "toggle switch: mutual repression" begin
            Congrego.extract_adjacency!(adj, [:pSrpR, :BM3R1, :pBM3R1, :SrpR]; pc=pc_qs2)
            @test adj[ni[:SrpR], ni[:BM3R1]] == Int8(-1)
            @test adj[ni[:BM3R1], ni[:SrpR]] == Int8(-1)
            @test adj[ni[:BM3R1], ni[:output]] == Int8(-1)
            @test count(!=(0), adj) == 3
        end

        @testset "duplicate TF: both inputs contribute" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :pLux, :SrpR]; pc=pc_qs2)
            @test adj[ni[:IPTG], ni[:SrpR]] == Int8(1)
            @test adj[ni[:AHL], ni[:SrpR]] == Int8(1)
            @test adj[ni[:IPTG], ni[:output]] == Int8(1)
            @test adj[ni[:AHL], ni[:output]] == Int8(1)
            @test count(!=(0), adj) == 4
        end

        # ── Group 7: Ordering invariance ─────────────────────────────
        @testset "OR/AND ordering collapse → same matrix" begin
            adj2 = zeros(Int8, N, N)
            Congrego.extract_adjacency!(adj, [:pTac, :TSrpR, :pBetI, :BM3R1]; pc=pc_qs2)
            Congrego.extract_adjacency!(adj2, [:pTac, :pBetI, :TSrpR, :BM3R1]; pc=pc_qs2)
            @test adj == adj2
        end

        @testset "debug printer + verbose hook run without error" begin
            Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :pSrpR]; pc=pc_qs2)
            @test (Congrego._print_adj_matrix(adj, pc_qs2); true)
            @test (Congrego.extract_adjacency!(adj, [:pTac, :SrpR, :pSrpR]; pc=pc_qs2, verbose=true); true)
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
