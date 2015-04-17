require 'gnuplot'

# Solves the 1D inear advection equation:
#
#   PDE: du/dt + a * du/dx = 0,
#   domain: 0 < x < xmax and 0 < t < tmax, 
#   periodic boundary consitions: u(0,t) = u(xmax, t),
#
# with stepsizes dx and dt.
#
# The initial condition must be supplied as a block, 
# and is saved as the proc object @ic.
class AdvectionEq

  attr_reader :dx, :dt, :a, :x, :t, :mx, :mt, :u, :method

  # Initialize the parameters and solve the equation
  # using one of the numeric solver schemes
  # :upwind, :lax_friedrichs, :leapfrog or :lax_wendroff
  # (default is :lax_wendroff).
  def initialize(params, &ic)
    params = { method: :lax_wendroff }.merge(params) 
    @dx = params['dx'] # x step size
    @dt = params['dt'] #t step size
    @a = params['a'] # parameter in the 1-dim linear advection eq.
    @ic = ic # initial condition
    @x = (0..params['xmax']).step(dx).to_a # x steps
    @t = (0..params['tmax']).step(dt).to_a # t steps
    @mx = @x.length # Number of x steps
    @mt = @t.length # Number of t steps
    @u = [] # Stores numerical solution

    @method = params[:method]
    case @method
      when :upwind then upwind
      when :lax_friedrichs then lax_friedrichs
      when :leapfrog then leapfrog
      when :lax_wendroff then lax_wendroff
      else raise(ArgumentError, "#{@method} is not a valid method.")
    end
  end

  # Evaluate the initial condition at x0
  def ic(x0)
    @ic.call(x0)
  end

  # Print the PDE
  def equation
    puts "PDE: du/dt + a * du/dx = 0, 0 < x < xmax and 0 < t < tmax, u(0,t) = u(xmax, t)"
  end

  # Plot the numeric solution
  def plot(title, *t_indices)
    return nil if u.empty?

    Gnuplot.open do |gp|
      Gnuplot::Plot.new(gp) do |plot|
        plot.title title.to_s
        plot.xlabel "x"
        plot.ylabel "u(x,t)"

        t_indices.each do |i|
          x = @x
          y = @u[i]
          plot.data << Gnuplot::DataSet.new([x,y]) do |ds|
            ds.with = "lines"
            ds.title = "t = #{@t[i]}"
          end
        end
      end
    end
  end

  ######### Available numeric schemes ##########
  
  private

    # Solves the 1D linear advection equation 
    # du/dt + a * du/dx = 0
    # with the Upwind scheme.
    # Requires: a > 0
    def upwind
      raise(ArgumentError, "the upwind scheme requires a>0.") if @a<=0
      @u[0] = @mx.times.map { |j| self.ic(@x[j]) } # IC: u(x,0)
      alpha = @a*@dt/@dx
      0.upto(@mt-2) do |n|
        @u[n+1] = (1.upto(@mx-2)).to_a
        @u[n+1].map! { |j| @u[n][j]-alpha*(@u[n][j]-@u[n][j-1])}
        @u[n+1].unshift @u[n][0]-alpha*(@u[n][0]-@u[n][@mx-2])
        @u[n+1].push @u[n+1][0] # periodic BC
      end
    end

    # Solves the 1D linear advection equation 
    # du/dt + a * du/dx = 0
    # with the Lachs-Friedrichs scheme.
    def lax_friedrichs
      @u[0] = @mx.times.map { |j| self.ic(@x[j]) } # IC: u(x,0)
      alpha = @a*@dt/@dx
      0.upto(@mt-2) do |n|
        @u[n+1] = (1.upto(@mx-3)).to_a
        @u[n+1].map! { |j| 0.5*(@u[n][j+1]+u[n][j-1])-0.5*alpha*(@u[n][j+1]-@u[n][j-1]) }
        @u[n+1].unshift 0.5*(@u[n][1]+u[n][@mx-2])-0.5*alpha*(@u[n][1]-@u[n][@mx-2]) # u(0,t)
        @u[n+1].push 0.5*(@u[n][0]+u[n][@mx-3])-0.5*alpha*(@u[n][0]-@u[n][@mx-3]) # u(max(x), t)
        @u[n+1].push @u[n+1][0] # periodic BC
      end
    end
    
    # Solves the 1D linear advection equation 
    # du/dt + a * du/dx = 0
    # with the Leapfrog scheme.
    def leapfrog
      @u[0] = @mx.times.map { |j| self.ic(@x[j]) } # IC: u(x,0)
      alpha = @a*@dt/@dx
      @u[1] = (1.upto(@mx-3)).to_a # u(x, t1)
      @u[1].map! { |j| 0.5*(@u[0][j+1]+u[0][j-1])-0.5*alpha*(@u[0][j+1]-@u[0][j-1]) }
      @u[1].unshift 0.5*(@u[0][1]+u[0][@mx-2])-0.5*alpha*(@u[0][1]-@u[0][@mx-2]) # u(0,t1)
      @u[1].push 0.5*(@u[0][0]+u[0][@mx-3])-0.5*alpha*(@u[0][0]-@u[0][@mx-3]) # u(max(x), t1)
      @u[1].push @u[1][0] # periodic BC
      1.upto(@mt-2) do |n|
        @u[n+1] = (1.upto(@mx-3)).to_a
        @u[n+1].map! { |j| @u[n-1][j]-alpha*(@u[n][j+1]-@u[n][j-1])}
        @u[n+1].unshift @u[n-1][0]-alpha*(@u[n][1]-@u[n][@mx-2]) # u(0,t)
        @u[n+1].push @u[n-1][@mx-2]-alpha*(@u[n][0]-@u[n][@mx-3]) # u(max(x), t)
        @u[n+1].push @u[n+1][0] # periodic BC
      end
    end

    # Solves the 1D linear advection equation 
    # du/dt + a * du/dx = 0
    # with the LaxWendroff scheme.
    def lax_wendroff
      @u[0] = @mx.times.map { |j| self.ic(@x[j]) } # IC: u(x,0)
      alpha = @a*@dt/@dx
      0.upto(@mt-2) do |n|
        @u[n+1] = (1.upto(@mx-3)).to_a
        @u[n+1].map! { |j| @u[n][j]-0.5*alpha*(@u[n][j+1]-@u[n][j-1])+0.5*(alpha**2)*(u[n][j+1]-2*u[n][j]+u[n][j-1]) }
        @u[n+1].unshift @u[n][0]-0.5*alpha*(@u[n][1]-@u[n][@mx-2])+0.5*(alpha**2)*(u[n][1]-2*u[n][0]+u[n][@mx-2]) # u(0,t)
        @u[n+1].push @u[n][@mx-2]-0.5*alpha*(@u[n][0]-@u[n][@mx-3])+0.5*(alpha**2)*(u[n][0]-2*u[n][@mx-2]+u[n][@mx-3]) # u(max(x), t)
        @u[n+1].push @u[n+1][0] # periodic BC
      end
    end
end
