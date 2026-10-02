#=
Adapted from part catalogue ms_roadblock_pTac.

The quorum sensing promoter :pLux and the quorum sensing signal :AHL have been added.
The input :aTC has been removed.

The pTra/TrarW parts have been removed because these require 3O8HSL induction

Removed transcription factors reported as toxic by (Nielsen et al 2016):
    :pQacR => :QacR,
    :pLitR => :LitR,
    :pAmtR => :AmtR,
    :pIcaRA => :IcaRA,

=#

module PC
using Congrego
export pc

const pc = PartCatalogue(
    parts=Set([
        :pTac,
        :pBAD,
        :pPhIF,
        :pBM3R1,
        :pPsrA,
        :pLmrA,
        :AraC,
        :PhIF,
        :PsrA,
        :BM3R1,
        :LmrA,
        :tfSp, # spacer
        :T10,  # strong terminator
        :pAmeR,
        :pHlyllR,
        # :pAmtR,
        :pBetI,
        # :pIcaRA,
        :pSrpR,
        # :pQacR,
        # :pLitR,
        :AmeR,
        :HlyllR,
        # :AmtR,
        # :IcaRA,
        :BetI,
        # :LitR,
        # :QacR,
        :SrpR,
        :TPhIF,  # inducible terminator
        :TBM3R1, # inducible terminator
        :TSrpR,  # inducible terminator
        :pLux,   # quorum sensing promoter
    ]),
    inputs=[:IPTG, :AHL],
    activatable_promoters=Dict(
        :pTac => :IPTG,
        :pLux => :AHL,
        :pBAD => :AraC,
    ),
    repressible_promoters=Dict(
        :pPhIF => :PhIF,
        :pBM3R1 => :BM3R1,
        :pPsrA => :PsrA,
        :pLmrA => :LmrA,
        :pAmeR => :AmeR,
        :pHlyllR => :HlyllR,
        # :pAmtR => :AmtR,
        # :pIcaRA => :IcaRA,
        :pBetI => :BetI,
        :pSrpR => :SrpR,
        # :pLitR => :LitR,
        # :pQacR => :QacR,
    ),
    terminators=Set([
        :T10,
    ]),
    roadblockers=Set([
        # activatable promoters
        :pTac,
        # repressible promoters
        :pPhIF,
        :pBM3R1,
        :pSrpR,
    ]),
    inducible_terminators=Dict{Symbol,Symbol}(
        :TPhIF => :PhIF,
        :TBM3R1 => :BM3R1,
        :TSrpR => :SrpR,
    ),
    part_length_map=Dict(
        :pTac => 135,
        :pLux => 125,
        :pBAD => 395,
        :AraC => 1056,
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
        # TODO - review these part lengths
        :TPhIF => 300,  # inducible terminator
        :TBM3R1 => 300, # inducible terminator
        :TSrpR => 467,  # inducible terminator
    ),
)

end # module PC