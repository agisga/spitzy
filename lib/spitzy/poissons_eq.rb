# Copyright (c) 2015 Alexej Gossmann 

# Numerically solves the 2D Poisson's equation:
#
#   * "Laplacian of u" = d^2 u / d^2 x + d^2 u / d^2 y = f(x,y),
#   * on the rectangular domain [xmin,xmax] x [ymin,ymax],
#   * with Dirichlet boundary conditions: u(x,y) = v(x,y) if 
#     (x,y) is on the boundary of the rectangular domain. 
#
class Poissons_eq

  # Attribute reader available for the following attributes:
  # * +x+      - array of the x coordinates of points in the domain 
  #              at which the numerical solution was evaluated
  # * +y+      - array of the y coordinates of points in the domain 
  #              at which the numerical solution was evaluated
  # * +mx+     - number of points on the side of the rectangular domain in x direction 
  # * +my+     - number of points on the side of the rectangular domain in y direction 
  # * +h+      - step size in x and y
  # * +mat+    - the matrix of the solved linear system
  # * +rhs+    - the right hand side of the solved linear system
  # * +method+ - the numerical scheme applied
  # * +u+      - the numerical solution as an array (the values are ordered as in +x+ and +y+). 
  #              It is obtained as a solution to a linear system.
  attr_reader :x, :y, :mx, :my, :h, :mat, :rhs, :method, :u

  # Constructor for all solver routines for the 2D Poisson's equation:
  #
  #   * "Laplacian of u" = d^2 u / d^2 x + d^2 u / d^2 y = f(x,y),
  #   * on the rectangular domain [xmin,xmax] x [ymin,ymax],
  #   * with Dirichlet boundary conditions: u(x,y) = v(x,y) if 
  #     (x,y) is on the boundary of the rectangular domain. 
  #   
  # It initializes the parameters and solves the equation
  # using one of the methods: 
  # :five_pt (5-point Laplacian), :nine_pt (9-point Laplacian)
  #
  # ==== Arguments
  #
  # * +xrange+  - An array of the form [xmin, xmax]. Defines the domain 
  #               on which the solution will be evaluated.
  #
  # * +yrange+  - An array of the form [ymin, ymax]. Defines the domain 
  #               on which the solution will be evaluated.
  #
  # * +h+       - The step size in x and y direction. It should be chosen such that
  #               subintervals of length +h+ cover the entire x range as well as the
  #               entire y range. 
  #
  # * +method+  - The numerical scheme used to solve the differential equation. Possible values are:
  #               +:five_pt+ (5-point Laplacian), +:nine_pt+ (9-point Laplacian)
  #               Default is +:five_pt+.
  #
  # * +bc+      - A +Proc+ or +Numeric+ object. The boundary condition v(x,y) of the differential 
  #               equation, which can be supplied as a +Proc+ object. If a +Numeric+ object is 
  #               supplied then v(x,y) is assumed to be constant equal to that number.
  #
  # * +f+       - A +Proc+ or +Numeric+ object. The right hand side f(x,y) of the differential 
  #               equation, which can be supplied as a +Proc+ object. If a +Numeric+ object is 
  #               supplied then f(x,y) is assumed to be constant equal to that number.
  #
  # ==== Usage
  #
  def initialize(xrange:, yrange:, mx:, my:, method: :five_pt, bc:, f:)
    raise(ArgumentError, "Expected xrange to be an array of length 2") unless xrange.length == 2
    raise(ArgumentError, "Expected yrange to be an array of length 2") unless yrange.length == 2
    @h = h
    @xmin, @xmax = xrange[0].to_f, xrange[1].to_f
    @ymin, @ymax = yrange[0].to_f, yrange[1].to_f

    xval = (@xmin..@xmax).step(@h).to_a 
    yval = (@ymin..@ymax).step(@h).to_a 
    @mx, @my = xval.length, yval.length
    unless xval[-1]==@xmax and yval[-1]==@ymax then
     raise(ArgumentError, "The step size +h+ should be chosen such that subintervals of 
                           length +h+ cover the entire x range as well as the entire y range.")
    end
    x, y = NMatrix::meshgrid([xval, yval], indexing: :ij)
    @x, @y = x.to_a, y.to_a
    # interior points
    @x_interior, @y_interior = NMatrix::meshgrid([xval[1..-2], yval[1..-2]], indexing: :ij)

    if f.is_a? Proc then
      @f = f 
    elsif f.is_a? Numeric then
      @f = Proc.new { |x| f }
    else
      raise(ArgumentError, "Expected Numeric or Proc input for a")
    end

    if bc.is_a? Proc then
      @bc = bc 
    elsif bc.is_a? Numeric then
      @bc = Proc.new { |x| bc }
    else
      raise(ArgumentError, "Expected Numeric or Proc input for a")
    end

    @method = method
    case @method
      when :five_pt then five_pt
      when :nine_pt then nine_pt
      else raise(ArgumentError, "#{@method} is not a valid method.")
    end
  end
  
  # Evaluate the function f(x,y) at x=x0 and y=y0
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  # * +y0+ - A floating point number
  #
  def f(x0, y0)
    @f.call(x0, y0)
  end

  # Evaluate the function bc(x,y) at x=x0 and y=y0
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  # * +y0+ - A floating point number
  #
  def bc(x0, y0)
    @bc.call(x0, y0)
  end


  ######### Available numeric schemes ##########
  
  private

    # Use the five-point Laplacian to solve the 2D Poisson's equation:
    #
    #   * "Laplacian of u" = d^2 u / d^2 x + d^2 u / d^2 y = f(x,y),
    #   * on the rectangular domain [xmin,xmax] x [ymin,ymax],
    #   * with Dirichlet boundary conditions: u(x,y) = v(x,y) if 
    #     (x,y) is on the boundary of the rectangular domain. 
    #   
    def five_pt
      # Number of interior points
      mx_interior = @my-2
      my_interior = @my-2
      m_interior = mx_interior*my_interior

      # Construct the right hand side for the linear system (only required at interiour points)
      @rhs = NMatrix.new([mx_interior, my_interior], dtype: :float64)
      @x_interior.each_with_indices do |v,i,j|
        if i==0 or i==mx_interior-1 or j==0 or j==my_interior-1 then
          rhsmat[i,j] = self.f(v, @y_interior[i,j]) - self.bc(v, @y_interior[i,j])
        else
          rhsmat[i,j] = self.f(v, @y_interior[i,j])
        end
      end
      @rhs.reshape!([m_interior, 1])

      # Construct the matrix for the linear system
      @mat = NMatrix.new([m_interior, m_interior], dtype: :float64)
      h2 = @h**2.0
      (0..m_interior).each {|i| @mat[i,i] = -4.0 / h2 }
      (0...m_interior).each {|i| @mat[i,i+1] = 1.0 / h2 }
      (0...m_interior).each {|i| @mat[i+1,i] = 1.0 / h2 }
      (0...(m_interior-my_interior)).each {|i| @mat[i,i+my_interior] = 1.0 / h2 }
      (0...(m_interior-my_interior)).each {|i| @mat[i+my_interior,i] = 1.0 / h2 }
      (1...mx_interior).each do |i|
        @mat[my_interior*i-1, my_interior*i] = 0
        @mat[my_interior*i, my_interior*i-1] = 0
      end

      # Compute the solution at the interior points
      u = @mat.solve(@rhs)
    end

end
