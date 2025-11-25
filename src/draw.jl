using Karnak
using Graphs
using Colors

function draw_pdf(ST::StratumTree, filename::String, width::Int64, height::Int64, cols::Vector{RGB{Float64}}, auts::Int64, draw_smoothing::Bool)
  nv = Graphs.nv(ST.GG)
  @pdf begin
      background("transparent")
      sethue("black")
      fontsize(8)
      drawgraph(ST.GG, 
          layout=stress, 
          vertexlabels = [ST.gen[v] == 0 ? "" : ST.gen[v] for v in 1:nv],
          vertexfillcolors = [cols[ST.col[v]] for v in 1:nv],
          margin = 13
      )
      info_string = string(auts)
      if draw_smoothing
        info_string = info_string * "," * string(length(ST.min_smoothings))
      end
      text(info_string, Point(-width/2, -height/2), halign=:left, valign = :top)
  end width height filename
end