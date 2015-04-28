# Copyright (c) 2015 Alexej Gossmann 

require 'nmatrix'

# Numerically solves a boundary value problem of the form:
#
#   * ODE: -a(x)*u''(x) + b(x)*u'(x) + c(x) = f(x), xmin<x<xmax,
#   * with boundary condition: u(xmin) = bc[1], u(xmax) = bc[2].
#   
class Bvp

  # Attribute reader available for the following attributes:
  # * +x+  - array of the points in the domain at which the numerical solution was evaluated
  # * +mx+ - number of points (i.e. length of +x+)
  # * +dx+ - step size in x
  # * +bc+ - the boundary conditions (an array of length two)
  # * +ax+ - value of a at every point in x
  # * +bx+ - value of b at every point in x
  # * +cx+ - value of c at every point in x
  # * +fx+ - value of f at every point in x
  # * +u+  - the numerical solution as an array
  # * +method+ - the numerical scheme applied
  attr_reader :x, :mx, :bc, :ax, :bx, :cx, :fx, :u, :method

  # Constructor for all solver routines for the boundary value problem:
  #
  #   * ODE: -a(x)*u''(x) + b(x)*u'(x) + c(x) = f(x), xmin<x<xmax,
  #   * with boundary condition: u(xmin) = bc[1], u(xmax) = bc[2].
  #
  # It initializes the parameters and solves the equation
  # using one of the methods: 
  # Linear Finite Elements. TODO: Implement more methods.
  #
  # ==== Arguments
  #
  # * +xrange+  - An array of the form [xmin, xmax]. The interval on which the solution will be evaluated.
  #
  # * +mx+      - Number of grid points. Determines the array x of equally spaced values on which the numerical
  #               solution can be evaluated.
  #
  # * +bc+      - An array of length two, specifying the two boundary conditions.
  #
  # * +method+  - The numerical scheme used to solve the BVP. Possible values are:
  #               +:lin_fin_elt+ (linear finite elements).
  #
  # * +&a+      - The function a(x) in the left hand side of the differential equation 
  #               which must be supplied as a +proc+ object. It is a function of x.
  #
  # * +&b+      - The function b(x) in the left hand side of the differential equation 
  #               which must be supplied as a +proc+ object. It is a function of x.
  #
  # * +&c+      - The function c(x) in the left hand side of the differential equation 
  #               which must be supplied as a +proc+ object. It is a function of x.
  #
  # * +&f+      - The right hand side of the differential equation which must be supplied as a +proc+ object.
  #               It is a function of x.
  #
  # ==== Usage
  #
  def initialize(xrange:, mx:, bc:, method: :lin_fin_elt, &a, &b, &c, &f)

    raise(ArgumentError, "Expected xrange to be an array of length 2") unless xrange.length == 2
    @mx = mx
    @xmin = xrange[0]
    @xmax = xrange[1]
    @dx = (@xmax-@xmin)/(@mx-1)
    @x = (@xmin..@xmax).step(@dx).to_a # x steps
    @bc = bc 
    @a = a 
    @b = b 
    @c = c 
    @f = f 
    @fa = [] # Stores the values of a(x) at every point in x
    @fb = [] # Stores the values of b(x) at every point in x
    @fc = [] # Stores the values of c(x) at every point in x
    @fx = [] # Stores the values of f(x) at every point in x
    @u = [] # Stores numerical solution

    @method = method
    case @method
      when :lin_fin_elt then lin_fin_elt
      else raise(ArgumentError, "#{@method} is not a valid method.")
    end
  end

  # Evaluate the function a(x) at x=x0 
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  def a(x0)
    @a.call(x0)
  end

  # Evaluate the function b(x) at x=x0
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  def b(x0)
    @b.call(x0)
  end

  # Evaluate the function c(x) at x=x0
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  def c(x0)
    @c.call(x0)
  end

  # Evaluate the function f(x) at x=x0
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  def f(x0)
    @f.call(x0)
  end

  ######### Available numeric schemes ##########
  
  private

    # Solve the initial value problem
    #  * dy/dx = f(x,y), xmin < x < xmax, 
    #  * y(xmin) = yini
    # with the forward Euler scheme u[n+1] = u[n] + dx * f[n].
    def lin_fin_el 
      asdf = NMatrix.new([2,2], [1,2,3,4])
      puts asdf
    end

end
