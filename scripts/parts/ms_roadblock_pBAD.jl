#=
Part catalogue used for Michiel's original investigation.
These produced a generator with 6 slots and 288 possible parts.

REVISIONS:
1) During the original investigation, the simulator was misconfigured
so that activatable promoters were not acting as roadblocks.
This behaviour has now been fixed in the simulator.
To replicate the behaviour during the original investigation,
the pBAD part was removed from the array of 'roadblockers' in ./ms_original.jl.


=#

module PC
using Congrego
export pc

const pc = PartCatalogue(
    parts=Set([
        :pTac,
        :pTet,
        :pBAD,
        :pTra,  # pTra* actually
        :pPhIF,
        :pBM3R1,
        :pPsrA,
        :pLmrA,
        :AraC,
        :TraRW,
        :PhIF,
        :PsrA,
        :BM3R1,
        :LmrA,
        :tfSp, # spacer
        :T10,  # strong terminator
        :pAmeR,
        :pHlyllR,
        :pAmtR,
        :pBetI,
        :pIcaRA,
        :pSrpR,
        :pQacR,
        :pLitR,
        :AmeR,
        :HlyllR,
        :AmtR,
        :IcaRA,
        :BetI,
        :LitR,
        :QacR,
        :SrpR,
        :TPhIF,  # inducible terminator
        :TBM3R1, # inducible terminator
        :TSrpR,  # inducible terminator
    ]),
    inputs=[:IPTG, :aTC],
    activatable_promoters=Dict(
        :pTac => :IPTG,
        :pTet => :aTC,
        :pBAD => :AraC,
        :pTra => :TraRW,
    ),
    repressible_promoters=Dict(
        :pPhIF => :PhIF,
        :pBM3R1 => :BM3R1,
        :pPsrA => :PsrA,
        :pLmrA => :LmrA,
        :pAmeR => :AmeR,
        :pHlyllR => :HlyllR,
        :pAmtR => :AmtR,
        :pIcaRA => :IcaRA,
        :pBetI => :BetI,
        :pSrpR => :SrpR,
        :pLitR => :LitR,
        :pQacR => :QacR,
    ),
    terminators=Set([
        :T10,
    ]),
    roadblockers=Set([
        :pBAD,
        :pPhIF,
        :pBM3R1,
        :pSrpR,
    ]),
    inducible_terminators=Dict(
        :TPhIF => :PhIF,
        :TBM3R1 => :BM3R1,
        :TSrpR => :SrpR,
    ),
    part_length_map=Dict(
        :pTac => 135,
        :pTet => 138,
        :pBAD => 395,
        :pTra => 119,
        :AraC => 1056,
        :TraRW => 847,
        :pAmeR => 127,
        :pHlyllR => 127,
        :pPhIF => 130,
        :pAmtR => 127,
        :pIcaRA => 131,
        :pPsrA => 130,
        :pBetI => 127,
        :pLitR => 127,
        :pQacR => 141,
        :pBM3R1 => 130,
        :pLmrA => 142,
        :pSrpR => 146,
        :AmeR => 838,
        :HlyllR => 788,
        :PhIF => 779,
        :AmtR => 847,
        :IcaRA => 743,
        :PsrA => 895,
        :BetI => 749,
        :LitR => 844,
        :QacR => 751,
        :BM3R1 => 751,
        :LmrA => 777,
        :SrpR => 823,
        :tfSp => 100,
        :T10 => 116,
        :TPhIF => 300,  # inducible terminator
        :TBM3R1 => 300, # inducible terminator
        :TSrpR => 300,  # inducible terminator
    ),
)

end # module PC