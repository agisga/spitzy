# Numerically solves an initial value problem of the form 
#
#  * dy/dx = f(x,y), xmin < x < xmax, 
#  * y(xmin) = yini
#   
class Ode

  # Attribute reader available for the following attributes:
  # * +x+  - array of the points in the domain at which the numerical solution was evaluated
  # * +mx+ - number of points (i.e. length of +x+)
  # * +fx+ - value of f at every point in x
  # * +u+  - the numerical solution as an array
  # * +tol+    - the error tolerance of the method (only applicable for methods
  #              with automatic step size adjustment)
  # * +method+ - the numerical scheme applied
  attr_reader :x, :mx, :fx, :u, :tol, :method

  # Constructor for all solver routines for the initial value problem:
  #  * dy/dx = f(x,y), xmin < x < xmax, 
  #  * y(xmin) = yini
  #
  # It initializes the parameters and solves the equation
  # using one of the methods: 
  # Forward Euler, Adams-Bashforth order 2, Adams-Moulton order 3, and more methods to be implemented
  #
  # ==== Arguments
  #
  # * +xrange+  - An array of the form [xmin, xmax]. The interval on which the solution will be evaluated.
  #
  # * +dx+      - Determines the grid of x values. The length of the x step for methods without automatic step size adjustment,
  #               or the maximum step size for methods with automatic step size adjustment. 
  #
  # * +yini+    - Initial value for y, i.e. y(xmin).
  #
  # * +tol+     - The desired error tolerance of the method (only for methods with automatic step size correction).
  # 
  # * +maxiter+ - The maximal number of performed iterations for methods with automatic step size adjustment.
  #
  # * +method+  - The numerical scheme used to solve the ODE. Possible values are:
  #               :dopri, :euler, :ab2. 
  #               Currently, the default is :dopri (Dormand-Prince is also currently the default method in MATLAB 
  #               and GNU Octave's ode45 solver and is the default choice for the Simulink's model explorer solver).
  #
  # * +&f+      - The right hand side of the differential equation which must be supplied as a +proc+ object.
  #               It is a function of x and y, where x should be the _first_ argument, and y the _second_.
  #
  # ==== Usage
  #
  # ...to be written
  #
  def initialize(xrange:, dx:, yini:, tol: 1e-2, maxiter: 1e6, method: :dopri, &f)

    raise(ArgumentError, "Expected xrange to be an array of length 2") unless xrange.length == 2
    @dx = dx
    @xmin = xrange[0]
    @xmax = xrange[1]
    @yini = yini
    @maxiter = maxiter
    @tol = tol
    @f = f 
    @x = [] # Stores the x grid
    @fx = [] # Stores the values of f at every point in x
    @u = [] # Stores numerical solution

    @method = method
    case @method
      when :euler then euler
      when :ab2 then ab2
      when :am3 then am3
      when :dopri then dopri
      else raise(ArgumentError, "#{@method} is not a valid method.")
    end
  end

  # Evaluate the function f(x,y) at x=x0 and y=y0 
  #
  # ==== Arguments
  #
  # * +x0+ - A floating point number
  #
  # * +y0+ - A floating point number representing y(x0)
  #
  def f(x0, y0)
    @f.call(x0, y0)
  end

  ######### Available numeric schemes ##########
  
  private

    # Solve the initial value problem
    #  * dy/dx = f(x,y), xmin < x < xmax, 
    #  * y(xmin) = yini
    # with the forward Euler scheme u[n+1] = u[n] + dx * f[n].
    def euler 
      @x = (@xmin..@xmax).step(@dx).to_a # x steps
      @x << @xmax if @x.last < @xmax
      @mx = @x.length # Number of x steps
      @u[0] = @yini
      @fx[0] = self.f(@x[0], @u[0])
      0.upto(@mx-2) do |n|
        @u[n+1] = @u[n] + @dx * @fx[n]
        @fx[n+1] = self.f(@x[n+1], @u[n+1])
      end
      #tol and maxiter not relevant for the Euler method
      @tol = nil
      @maxiter = nil
    end
   
    # Solve the initial value problem
    #  * dy/dx = f(x,y), xmin < x < xmax, 
    #  * y(xmin) = yini
    # with the Adams-Bashforth method of order 2, 
    # given by the explicit formula u[n+1] = u[n] + dx/2 * (3*f[n] - f[n-1]).
    # Since the method requires the last two functional values for the approximation 
    # of the next functional value, a Runge-Kutta method (Heun’s method) of 
    # order 2 is applied to approximate the functional value at the second time step.
    #
    def ab2 
      @x = (@xmin..@xmax).step(@dx).to_a # x steps
      @x << @xmax if @x.last < @xmax
      @mx = @x.length # Number of x steps

      @u[0] = @yini
      @fx[0] = self.f(@x[0], @u[0])
      # Runge-Kutta of order 2 for the second time step
      @u[1] = @u[0] + @dx/2.0 * (@fx[0] + self.f(@x[1], @u[0] + @dx*@fx[0]))
      @fx[1] = self.f(@x[1], @u[1])

      1.upto(@mx-2) do |n|
        @u[n+1] = @u[n] + @dx/2.0 * (3.0*@fx[n] - @fx[n-1])
        @fx[n+1] = self.f(@x[n+1], @u[n+1])
      end

      #tol and maxiter not relevant for the this method
      @tol = nil
      @maxiter = nil
    end
 
    # Solve the initial value problem
    #  * dy/dx = f(x,y), xmin < x < xmax, 
    #  * y(xmin) = yini
    # with the forward Dormand-Prince method.
    #
    # This method automatically adapts the step size in order to keep the error
    # of the numerical solution below the tolerance level +tol+. However, the step size
    # is not allowed to exceed the specified maximal step size +dx+.
    # The algorithm throws an exception if it fails to compute the numerical solution
    # on the given x domain in less than or equal to +maxiter+ iterations.
    # 
    # === Reference
    #
    #     J.R. Dormand, P.J. Prince, "A family of embedded Runge-Kutta formulae"
    #
    #     E.Hairer, S.P.Norsett and G.Wanner, "Solving Differential Equations I, Nonstiff Problems", p. 176.
    # 
    def dopri
      a21 = 1.0/5.0
      a31 = 3.0/40.0
      a32 = 9.0/40.0
      a41 = 44.0/45.0
      a42 = -56.0/15.0
      a43 = 32.0/9.0
      a51 = 19372.0/6561.0
      a52 = -25360.0/2187.0
      a53 = 64448.0/6561.0
      a54 = -212.0/729.0
      a61 = 9017.0/3168.0
      a62 = -355.0/33.0
      a63 = 46732.0/5247.0
      a64 = 49.0/176.0
      a65 = -5103.0/18656.0
      a71 = 35.0/384.0
      a72 = 0.0
      a73 = 500.0/1113.0
      a74 = 125.0/192.0
      a75 = -2187.0/6784.0
      a76 = 11.0/84.0

      c2 = 1.0 / 5.0
      c3 = 3.0 / 10.0
      c4 = 4.0 / 5.0
      c5 = 8.0 / 9.0
      c6 = 1.0
      c7 = 1.0

      b1order5 = 35.0/384.0
      b2order5 = 0.0
      b3order5 = 500.0/1113.0
      b4order5 = 125.0/192.0
      b5order5 = -2187.0/6784.0
      b6order5 = 11.0/84.0
      b7order5 = 0.0

      b1order4 = 5179.0/57600.0
      b2order4 = 0.0
      b3order4 = 7571.0/16695.0
      b4order4 = 393.0/640.0
      b5order4 = -92097.0/339200.0
      b6order4 = 187.0/2100.0
      b7order4 = 1.0/40.0

      @x[0] = @xmin
      @u[0] = @yini
      @fx[0] = self.f(@x[0], @u[0])
      h = @dx 
      i = 0

      0.upto(@maxiter) do |iter|
         # Compute the function values
         k1 = @fx[i] 
         k2 = self.f(@x[i] + c2*h, @u[i] + h*(a21*k1))
         k3 = self.f(@x[i] + c3*h, @u[i] + h*(a31*k1+a32*k2))
         k4 = self.f(@x[i] + c4*h, @u[i] + h*(a41*k1+a42*k2+a43*k3))
         k5 = self.f(@x[i] + c5*h, @u[i] + h*(a51*k1+a52*k2+a53*k3+a54*k4))
         k6 = self.f(@x[i] + h, @u[i] + h*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
         k7 = self.f(@x[i] + h, @u[i] + h*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6))

         error = (b1order5 - b1order4)*k1 + (b3order5 - b3order4)*k3 + (b4order5 - b4order4)*k4 + 
           (b5order5 - b5order4)*k5 + (b6order5 - b6order4)*k6 + (b7order5 - b7order4)*k7
         error = error.abs

         # error control
         if error < @tol then
            @x[i+1] = @x[i] + h
            @u[i+1] = @u[i] + h * (b1order5*k1 + b3order5*k3 + b4order5*k4 + b5order5*k5 + b6order5*k6)
            @fx[i+1] = self.f(@x[i+1], @u[i+1])
            i = i+1
         end

         delta = 0.84 * (@tol / error)**0.2
         if delta <= 0.1 then
            h = h * 0.1
         elsif delta >= 4.0 then
            h = h * 4.0
         else 
            h = delta * h
         end

         # set h to the user specified maximal allowed value
         h = @dx if h > @dx 

         if @x[i] >= @xmax then
           break
         elsif @x[i] + h > @xmax then
           h = @xmax - @x[i]
         end
      end

      @mx = @x.length # Number of x steps

      raise(RuntimeError, "Maximal number of iterations reached 
            before evaluation of the solution on the entire x interval 
            was completed (try to increase maxiter or use a different method") if @x.last < @xmax
    end
 
end
