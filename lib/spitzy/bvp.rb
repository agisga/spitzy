# Copyright (c) 2015 Alexej Gossmann 

# Numerically solves a boundary value problem of the general form:
#
#   * ODE: -(a*u')' + b*u' + c*u = f, xmin<x<xmax,
#   * where a, b, c, f and u are functions of x,
#   * with Dirichlet boundary conditions: u(xmin) = bc[1], u(xmax) = bc[2].
#   
class Bvp

  # Attribute reader available for the following attributes:
  # * +x+      - array of the points in the domain at which the numerical solution was evaluated
  # * +mx+     - number of points (i.e. length of +x+)
  # * +dx+     - step size in x
  # * +bc+     - the boundary conditions (an array of length two)
  # * +method+ - the numerical scheme applied
  # * +u+      - the numerical solution as an array. It is obtained as a solution to a linear system.
  # * +rhs+    - the right hand side of the solved linear system
  # * +mat+    - the stiffness matrix, i.e. the matrix of the solved linear system
  attr_reader :x, :mx, :bc, :method, :u, :rhs, :mat

  # Constructor for all solver routines for the boundary value problem:
  #
  #   * ODE: -(a*u')' + b*u' + c*u = f, xmin<x<xmax,
  #   * where a, b, c, f and u are functions of x,
  #   * with Dirichlet boundary conditions: u(xmin) = bc[1], u(xmax) = bc[2].
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
  # * +a+       - A +Proc+ or +Numeric+ object. The function a(x) in the left hand side 
  #               of the differential equation which can be supplied as a +Proc+ object.
  #               If a +Numeric+ object is supplied then a(x) is assumed to be constant equal to that number.
  #
  # * +b+       - A +Proc+ or +Numeric+ object. The function b(x) in the left hand side 
  #               of the differential equation which can be supplied as a +Proc+ object.
  #               If a +Numeric+ object is supplied then b(x) is assumed to be constant equal to that number.
  #
  # * +c+       - A +Proc+ or +Numeric+ object. The function c(x) in the left hand side 
  #               of the differential equation which can be supplied as a +Proc+ object.
  #               If a +Numeric+ object is supplied then c(x) is assumed to be constant equal to that number.
  #
  # * +f+       - A +Proc+ or +Numeric+ object. The right hand side f(x) of the differential 
  #               equation, which can be supplied as a +Proc+ object. If a +Numeric+ object is 
  #               supplied then f(x) is assumed to be constant equal to that number.
  #
  # ==== Usage
  #
  def initialize(xrange:, mx:, bc:, method: :lin_fin_elt, a:, b:, c:, f:)

    raise(ArgumentError, "Expected xrange to be an array of length 2") unless xrange.length == 2
    raise(ArgumentError, "Expected bc to be an array of length 2") unless bc.length == 2
    @mx = mx
    @xmin = xrange[0].to_f
    @xmax = xrange[1].to_f
    @dx = (@xmax-@xmin) / (@mx-1)
    @x = (@xmin..@xmax).step(@dx).to_a # x steps
    @bc = bc 
    @u = [] # Stores numerical solution

    if a.is_a? Proc then
      @a = a 
    elsif a.is_a? Numeric then
      @a = Proc.new { |x| a }
    else
      raise(ArgumentError, "Expected Numeric or Proc input for a")
    end
    if b.is_a? Proc then
      @b = b 
    elsif b.is_a? Numeric then
      @b = Proc.new { |x| b }
    else
      raise(ArgumentError, "Expected Numeric or Proc input for b")
    end
    if c.is_a? Proc then
      @c = c 
    elsif c.is_a? Numeric then
      @c = Proc.new { |x| c }
    else
      raise(ArgumentError, "Expected Numeric or Proc input for c")
    end
    if f.is_a? Proc then
      @f = f 
    elsif f.is_a? Numeric then
      @f = Proc.new { |x| f }
    else
      raise(ArgumentError, "Expected Numeric or Proc input for f")
    end

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

    # Solve the boundary value problem
    #
    #   * ODE: -(a*u')' + b*u' + c*u = f, xmin<x<xmax,
    #   * where a, b, c, f and u are functions of x,
    #   * with Dirichlet boundary conditions: u(xmin) = bc[1], u(xmax) = bc[2]
    #
    # using the finite elements method with piecewise linear elements.
    #
    # ==== References
    #
    # * A. Quarteroni, R. Sacco, F. Saleri, "Numerical Mathematics", 2nd ed., section 7.4.
    #
    def lin_fin_elt 
      gridcenters = ((@xmin+@dx/2)..(@xmax-@dx/2)).step(@dx).to_a
      a_gridcenters = NMatrix.new([@mx-1,1], gridcenters.map {|x0| self.a(x0) }, dtype: :float64)
      b_gridcenters = NMatrix.new([@mx-1,1], gridcenters.map {|x0| self.b(x0) }, dtype: :float64) 
      c_gridcenters = NMatrix.new([@mx-1,1], gridcenters.map {|x0| self.c(x0) }, dtype: :float64) 
      f_gridcenters = NMatrix.new([@mx-1,1], gridcenters.map {|x0| self.f(x0) }, dtype: :float64) 

      # Construct the  RHS for the linear system
      @rhs = (f_gridcenters[0..-2] + f_gridcenters[1..-1]) * @dx * 0.5
      @rhs[0] = @rhs[0] - @bc[0] * (-a_gridcenters[0] / @dx - b_gridcenters[0] / 2.0 + @dx * c_gridcenters[0] / 3.0)
      @rhs[-1] = @rhs[-1] - @bc[1] * (-a_gridcenters[-1] / @dx + b_gridcenters[-1] / 2.0 + @dx * c_gridcenters[-1] / 3.0)
      
      # Construct the tridiagonal stiffness matrix
      # TODO: the matrix should be saved as a sparse stype
      dd = NMatrix.new([@mx-2,1], 0.0, dtype: :float64) 
      dc = NMatrix.new([@mx-2,1], 0.0, dtype: :float64) 
      dr = NMatrix.new([@mx-2,1], 0.0, dtype: :float64)  
      ld = NMatrix.new([@mx-2,1], 0.0, dtype: :float64)  
      lc = NMatrix.new([@mx-2,1], 0.0, dtype: :float64)  
      lr = NMatrix.new([@mx-2,1], 0.0, dtype: :float64)  
      ud = NMatrix.new([@mx-2,1], 0.0, dtype: :float64)  
      uc = NMatrix.new([@mx-2,1], 0.0, dtype: :float64)  
      ur = NMatrix.new([@mx-2,1], 0.0, dtype: :float64) 
      1.upto(@mx-2) do |i|
        dd[i-1] = (a_gridcenters[i-1] + a_gridcenters[i]) / @dx
        dc[i-1] = (b_gridcenters[i-1] - b_gridcenters[i]) / 2.0
        dr[i-1] = @dx * (c_gridcenters[i-1] + c_gridcenters[i]) / 3.0
        if i>1 then
          ld[i-2] = -a_gridcenters[i-1] / @dx
          lc[i-2] = -b_gridcenters[i-1] / 2.0
          lr[i-2] = @dx * c_gridcenters[i-1] / 6.0
        end
        if i<(@mx-2) then
          ud[i-1] = -a_gridcenters[i] / @dx
          uc[i-1] = b_gridcenters[i] / 2.0
          ur[i-1] = @dx * c_gridcenters[i] / 6.0
        end
      end
      diag = dd + dc + dr
      ldiag = ld + lc + lr
      udiag = ud + uc + ur
      @mat = NMatrix.diagonal(diag.to_a, dtype: :float64)
      1.upto(@mx-3) do |j|
        @mat[j,j-1] = ldiag[j-1]
        @mat[j-1,j] = udiag[j-1]
      end

      # Solve the system
      # TODO: The system is tridiagonal and should be solved accordingly
      w = @mat.solve(@rhs)
      @u = w.to_a
      @u.insert(0, @bc[0])
      @u.insert(-1, @bc[1])
    end

end
