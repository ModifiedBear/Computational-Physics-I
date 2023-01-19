using LinearAlgebra
using CairoMakie
using BenchmarkTools

function get_boundary_indices(nx, ny, bsize)
  # this can be improved 
  idx = falses(nx, ny)
  idx[:,1:bsize] .= true
  idx[:,ny-bsize+1:ny] .= true
  idx[1:bsize,:] .= true
  idx[nx-bsize+1:nx,:] .= true

  #return idx
  #return findall(vec(idx)) # return linear indices
  return getindex.(findall(idx), [[1] [2]])
end

function init(nx, ny)
  # initial values
  L = 10  
  #  calculations
  #dx, dy = L/nx, L/ny
  xs = LinRange(-L,L,nx)
  ys = LinRange(-L,L,ny)
  
  
  
      
  c = 0.3
  h = 1. # spatial width
  k = 1. # time step width
  α = ones(nx, ny) .* ((c*k) / h)^2 # alpha squared ()
  κ = α .* (k/h);

  #α[134:138, 1:45] .= 0;
  #α[130:134,:] .= 0;
  #α[130:134,70:80] .= ((c*k) / h)^2;
  #α[134:138, 1:45] .= 0;
  #α[134:138, 55:95] .= 0;
  #α[134:138, 105:end] .= 0;
  
  Z_new = [exp(- 0.1 * (x-5)^2 - 0.1 * (y-5)^2) for x in xs, y in ys]
  Z = zeros(nx, ny)
  Z_old = zeros(nx, ny)
  #data = zeros(N, nx+1, ny+1)
  
  
  return Z, Z_new, Z_old, α, κ, xs, ys
end

function oh2(u, u_new, u_old, α, κ, M, nx, ny,n)
  # better than the vectorized format
  boundary=get_boundary_indices(nx, ny, 1)
  for ii in 1:n
      # you NEED to use copy(), this isn't python
      ω = 0.15
      #u_new[nx-2:nx-1,:] .= 10 .* sin(ω * ii)
      M[ii,:,:] = copy(u_new)
      u_old = copy(u)
      u = copy(u_new)
      for ii in 2:nx-1
          for jj in 2:ny-1
              u_new[ii, jj]  = α[ii,jj] * (u[ii-1, jj] + u[ii+1, jj] + u[ii, jj-1] + u[ii, jj+1] - 4*u[ii, jj])                 
              u_new[ii, jj] += 2 * u[ii, jj] - u_old[ii, jj]

              # absorbing boundaries
              #u_new[1, jj] = u[2,   jj]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (u_new[2,   jj]-u[1, jj])
              #u_new[nx,jj] = u[nx-1,jj]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (u_new[nx-1,jj]-u[nx,jj])
              #u_new[ii, 1] = u[ii,   2]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (u_new[ii,   2]-u[ii, 1])
              #u_new[ii,ny] = u[ii,ny-1]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (u_new[ii,ny-1]-u[ii,ny])
          end
      end
      # no need to update boudnary
      #u_new[boundary] .= 0
  end
  
  return M
end  

function get_plots(V, x, y)
  set_theme!(theme_black())
  custom_cmap=[RGBAf(c.r, c.g, c.b, 0.05) for c in to_colormap(:grays)]
  fig = Figure(resolution=(500,500))
  ax = Axis(fig[1,1])
  ax.aspect = DataAspect() 
  heatmap!(ax,x,y, norm.(V); colormap=Reverse(:Spectral), interpolate=false)
  fig
end

begin
  N = 200
  nx, ny = 100,100
  U, U_new, U_old, A,K, xs, ys = init(nx,ny);
  data = oh2(U, U_new, U_old, A, K, zeros(N, nx, ny),nx, ny, N);
  get_plots(data[end,:,:], xs, ys)
  
end

begin
  set_theme!(theme_black())
  custom_cmap=[RGBAf(c.r, c.g, c.b, 0.05) for c in to_colormap(:grays)]
  f = Figure(resolution=(500,500))   
  ax = Axis(f[1,1])
  ax.aspect=DataAspect()
  rowsize!(f.layout, 1, ax.scene.px_area[].widths[2]) # set colorbar height
  framerate = 30
  record(f, "wave_2_non_absorbing.mp4", 1:2:N;
          framerate = framerate) do idx
      hm = heatmap!(ax, xs, ys, norm.(data[idx,:,:]); colormap=Reverse(:Spectral), colorrange=(0,maximum(data)), interpolate=true)
      
      #bl = heatmap!(ax,alpha[diagind(alpha)],colormap = :grays,colorrange=(0,0), transparency=true)
      bl = heatmap!(ax, xs, ys, A; colormap=custom_cmap, interpolate=true)
      
  end
end
