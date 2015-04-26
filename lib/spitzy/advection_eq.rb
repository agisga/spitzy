# Copyright (c) 2015 Alexej Gossmann 


# Numerically solves the 1D linear advection equation:
#
#  * PDE: du/dt + a * du/dx = 0,
#  * on the domain: 0 < x < xmax and 0 < t < tmax, 
#  * with periodic boundary consitions: u(0,t) = u(xmax, t),
#  * with initial condition as supplied by the user.
#   
class AdvectionEq

  # Attribute reader available for the following attributes:
  # * +dx+ - time step in x                   
  # * +dt+ - time step in t
  # * +a+  - parameter in the 1D linear advection equation
  # * +x+  - array of space points
  # * +t+  - array of time points
  # * +mx+ - number of space points (i.e. length of +x+)
  # * +mt+ - number of time points (i.e. length of +t+)
  # * +u+  - the numerical solution as an array of arrays
  # * +method+ - the numerical scheme applied
  attr_reader :dx, :dt, :a, :x, :t, :mx, :mt, :u, :method

  # Constructor for all solver routines for the 1D linear advection equation:
  #  * PDE: du/dt + a * du/dx = 0,
  #  * on the domain: xmin < x < xmax and tmin < t < tmax, 
  #  * with periodic boundary consitions: u(xmin,t) = u(xmax,t),
  #  * with initial condition u(x,tmin) as supplied by the user.
  #
  # It initializes the parameters and solves the equation
  # using one of the methods: 
  # Upwind, Leapfrog, Lax-Friedrichs, Lax-Wendroff.
  #
  # ==== Arguments
  #
  # * +xrange+  - An array of the form [xmin, xmax]. The space domain xmin<=x<=xmax 
  #               on which the solution will be evaluated.
  #
  # * +trange+  - An array of the form [tmin, tmax]. The time domain tmin<=t<=tmax 
  #               on which the solution will be evaluated.
  #
  # * +dx+      - Stepsize in space, i.e. in x
  #
  # * +dt+      - Stepsize in time, i.e. in t
  #
  # * +a+       - The constant speed +a+ in the PDE du/dt + a * du/dx = 0
  #
  # * +method+ - The numerical scheme used to solve the PDE. Possible values are:
  #              :upwind, :lax_friedrichs, :leapfrog or :lax_wendroff (default is :lax_wendroff).
  #
  # * +&ic+     - The initial condition which must be supplied as a +proc+ object.
  #               The initial condition is a function in x at time t=0.
  #
  # ==== Usage
  #
  # ic = proc { |x| Math::cos(2*Math::PI*x) + 0.2*Math::cos(10*Math::PI*x) }
  # numsol = AdvectionEq.new(xrange: [0.0,1.0], trange: [0.0,10.0], dx: 1.0/101, 
  #                          dt: 0.95/101, a: 1.0, method: :upwind, &ic)
  #
  def initialize(xrange: , trange: , dx: , dt: , a: , method: :lax_wendroff, &ic)
    raise(ArgumentError, "Expected xrange to be an array of length 2") unless xrange.length == 2
    raise(ArgumentError, "Expected trange to be an array of length 2") unless trange.length == 2
    @dx = dx # x step size
    @dt = dt #t step size
    @a = a # parameter in the 1-dim linear advection eq.
    @ic = ic # initial condition
    @x = (xrange[0]..xrange[1]).step(dx).to_a # x steps
    @x << xrange[1] if @x.last < xrange[1]
    @mx = @x.length # Number of x steps
    @t = (trange[0]..trange[1]).step(dt).to_a # t steps
    @t << trange[1] if @t.last < trange[1]
    @mt = @t.length # Number of t steps
    @u = [] # Stores numerical solution

    @method = method
    case @method
      when :upwind then upwind
      when :lax_friedrichs then lax_friedrichs
      when :leapfrog then leapfrog
      when :lax_wendroff then lax_wendroff
      else raise(ArgumentError, "#{@method} is not a valid method.")
    end
  end

  # Evaluate the initial condition at x0
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  def ic(x0)
    @ic.call(x0)
  end

  # Get a character string representing the PDE, its domain and boundary condition. 
  def equation
    "du/dt + #{@a} * du/dx = 0, #{@x[0]} < x < #{@x[-1]} and #{@t[0]} < t < #{@t[-1]}, u(0,t) = u(#{@x[-1]}, t)"
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
    # with the Lax-Friedrichs scheme.
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
    # with the Lax-Wendroff scheme.
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
