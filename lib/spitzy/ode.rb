# Numerically solves an initial value problem of the form 
#
#  * dy/dt = f(t,y), tmin < t < tmax, 
#  * y(tmin) = yini
#   
class Ivp

  # Attribute reader available for the following attributes:
  # * +t+  - array of the points in the domain at which the numerical solution was evaluated
  # * +mt+ - number of points (i.e. length of +t+)
  # * +ft+ - value of f at every point in t
  # * +u+  - the numerical solution as an array
  # * +method+ - the numerical scheme applied
  attr_reader :t, :mt, :ft, :u, :method

  # Constructor for all solver routines for the initial value problem:
  #  * dy/dt = f(t,y), tmin < t < tmax, 
  #  * y(tmin) = yini
  #
  # It initializes the parameters and solves the equation
  # using one of the methods: 
  # Forward Euler, Adams-Bashforth order 2, Adams-Moulton order 3, and more methods to be implemented
  #
  # ==== Arguments
  #
  # * +t+       - An array of points at which the numerical solution will be evaluated.
  #               Ideally, the points should be equally spaced.
  #
  # * +yini+    - Initial value for y, i.e. y(t[0])
  #
  # * +:method+ - The numerical scheme used to solve the ODE. Possible values are:
  #               :euler, :ab2, :am3 (default is :am3).
  #
  # * +&f+      - The right hand side of the differential equation which must be supplied as a +proc+ object.
  #               It is a function of y and t, where y should be the _first_ argument, and t the _second_.
  #
  # ==== Usage
  #
  # ...to be written
  #
  def initialize(t:, yini:, method:, &f)
    @yini = yini
    @t = t
    @f = f 

    @mt = @t.length # Number of t steps
    @ft = [] # Stores the values of f at every point in t
    @u = [] # Stores numerical solution

    @method = method
    case @method
      when :euler then euler
      when :ab2 then ab2
      when :am3 then am3
      else raise(ArgumentError, "#{@method} is not a valid method.")
    end
  end

  # Evaluate the function f(y,t) at y=y0 and t=t0
  #
  # ==== Arguments
  #
  # * +t0+ - A floating point number
  #
  # * +y0+ - A floating point number representing y(t0)
  #
  def f(y0, t0)
    @f.call(y0, t0)
  end

  ######### Available numeric schemes ##########
  
  private

    # Solve the initial value problem
    #  * dy/dt = f(t,y), tmin < t < tmax, 
    #  * y(tmin) = yini
    # with the forward Euler scheme.
    def euler 
      @u[0] = @yini
      @ft[0] = self.f(@u[0], @t[0])
      0.upto(@mt-2) do |n|
        @u[n+1] = @u[n] + (@t[n+1] - @t[n]) * @ft[n]
        @ft[n+1] = self.f(@u[n+1], @t[n+1])
      end
    end
end
