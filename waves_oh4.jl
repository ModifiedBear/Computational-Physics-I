using LinearAlgebra
using CairoMakie
using BenchmarkTools

function init(nx, ny)
  # initial values
  L = 10  
  #  calculations
  #dx, dy = L/nx, L/ny
  xs = LinRange(-L,L,nx);
  ys = LinRange(-L,L,ny);
  
  
  
      
  c = 0.3;
  h = 0.5 ;# spatial width
  k = 0.2; # time step width (remember, dt < dx^2/2)?
  α = ones(nx, ny) .* ((c*k) / h)^2; # alpha squared ()
  #α[144:148,  1:55] .= 0;
  #α[144:148, 60:90] .= 0;
  #α[144:148, 95:ny] .= 0;

  #α[div(nx,3):nx-div(nx,3),div(ny,3):ny-div(ny,3)] .= 0


  κ = sqrt.(α) .* (k/h) .* 10; # no tengo idea por qué hay que multiplicar por 10 para que no refleje

 
  Zₙ₊₁ = zeros(nx, ny);#[exp(- 0.1 * (x-5)^2 - 0.1 * (y-5)^2) for x in xs, y in ys]
  Zₙ   = zeros(nx, ny);
  Zₙ₋₁ = zeros(nx, ny);
  #data = zeros(N, nx+1, ny+1)
  
  
  return Zₙ, Zₙ₊₁, Zₙ₋₁, α, κ, xs, ys;
end

function oh4(vₙ, vₙ₊₁, α, κ, nx, ny,n)
  # better than the vectorized format
  mat = zeros(n, nx, ny)
  
  bndry = 2
  for ii in 1:n
    mat[ii,:,:] = vₙ
      # you NEED to use copy(), this isn't python
    ω = 2* π * 0.02;
    #vₙ₊₁[nx-bndry-1:nx-bndry,div(ny,3):ny-div(ny,3)] .= 10 .* cos(ω * ii);
    #vₙ₊₁[nx-bndry-2:nx-bndry,ny-bndry-2:ny-bndry] .= 10 .* cos(ω * ii);
    vₙ₊₁[div(nx,2)-1:div(nx,2)+1,div(ny,2)-1:div(ny,2)+1] .= 10 .* cos(ω * ii);
    #M[ii,:,:] = copy(vₙ₊₁);
    vₙ₋₁ = copy(vₙ);
    vₙ = copy(vₙ₊₁);
    
    for ii in 3:nx-2
      for jj in 3:ny-2
        vₙ₊₁[ii, jj]  = α[ii,jj] * (  -vₙ[ii,jj-2]+
                                    16*vₙ[ii,jj-1]-
                                       vₙ[ii-2,jj]+
                                    16*vₙ[ii-1,jj]-
                                    60*vₙ[ii,  jj]+
                                    16*vₙ[ii+1,jj]-
                                       vₙ[ii+2,jj]+
                                    16*vₙ[ii,jj+1]-
                                       vₙ[ii,jj+2]
                                    ) + 
                            2 * vₙ[ii, jj] - vₙ₋₁[ii, jj];
        
          # absorbing boundaries
          # don't mix for loops and vectorization, bad results
        for kk in 1:bndry
          vₙ₊₁[kk, jj] =     vₙ[kk+1,   jj]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[kk+1,   jj]-vₙ[kk,     jj]);# x = 0
          vₙ₊₁[nx-2+kk,jj] = vₙ[nx-3+kk,jj]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[nx-3+kk,jj]-vₙ[nx-2+kk,jj]);# x = N
          vₙ₊₁[ii, kk] =     vₙ[ii,   kk+1]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[ii,   kk+1]-vₙ[ii,     kk]);# y = 0
          vₙ₊₁[ii,ny-2+kk] = vₙ[ii,ny-3+kk]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[ii,ny-3+kk]-vₙ[ii,ny-2+kk]);# y = N
        end
      end
    end
  end

  return mat;
end

function get_plot(data,medium, xs, ys)
  custom_cmap=[RGBAf(0, 0, 0, 1) for c in to_colormap(:Greys)]
  fig=Figure()
  ax=Axis(fig[1,1])
  ax.aspect = DataAspect() 
  hm=Makie.heatmap!(ax,xs, ys, data,colormap=:Spectral, interpolate=false)
  
  mask = ones(size(medium))
  mask[medium .> minimum(medium)] .= NaN;
  
  #bl=Makie.heatmap!(ax, xs, ys, mask; colormap=custom_cmap, interpolate=false)
  Colorbar(fig[1,2],hm)
  rowsize!(fig.layout, 1, ax.scene.px_area[].widths[2]) # set colorbar height
  fig
end

begin
  N = 330;
  nx, ny = 151,151;
  U, U_new, U_old, Alpha,Kappa, xs, ys = init(nx,ny);
  #@benchmark data = oh4(U, U_new, U_old, A, K, zeros(N, nx, ny),nx, ny, N)
  data = oh4(U, U_new, Alpha, Kappa, nx, ny, N);
  get_plot(norm.(data[end,:,:]),Alpha,xs,ys)
end


