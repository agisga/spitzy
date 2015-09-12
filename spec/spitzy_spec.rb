require 'spec_helper'

describe Spitzy do
  it 'has a version number' do
    expect(Spitzy::VERSION).not_to be nil
  end

  describe Spitzy::AdvectionEq do

    #  solve the 1D linear advection equation given as:
    #  * PDE: $\frac{du}{dt} + a \frac{du}{dx} = 0$,
    #  * on the domain: $0 < x < 1$ and $0 < t < 10$, 
    #  * with periodic boundary consitions: $u(0,t) = u(1, t)$,
    #  * with initial condition: $u(x,0) = \cos(2\pi x) + \frac{1}{5}\cos(10\pi x)$.

    [:upwind, :lax_friedrichs, :leapfrog, :lax_wendroff].each do |mthd|
      context "when method: #{mthd}" do
        subject(:numsol) do
          ic = proc { |x| Math::cos(2*Math::PI*x) + 0.2*Math::cos(10*Math::PI*x) }

          Spitzy::AdvectionEq.new(xrange: [0.0,1.0], trange: [0.0, 10.0], 
                                  dx: 1.0/1001, dt: 0.95/1001, a: 1.0,
                                  method: mthd, &ic)
        end

        it "reports the utilized method" do
          expect(numsol.method).to eq(mthd)
        end

        it "computes a precise numerical solution" do
          exactsol = proc { |x,t| Math::cos(2*Math::PI*(x-t)) + 0.2*Math::cos(10*Math::PI*(x-t)) }
          u_exact = []
          numsol.t.each { |t| u_exact << (0...numsol.mx).each.map { |i| exactsol.call(numsol.x[i],t) } }
          u_exact.flatten!
          error = []
          u_num = numsol.u.flatten
          (0...(numsol.mx*numsol.mt)).each { |i| error << (u_num[i] - u_exact[i]).abs }
          maxerror = error.max

          expect(maxerror).to be_within(0.1).of(0.0)
        end
      end 
    end
  end

  describe Spitzy::Bvp do

    #  solve the following boundary value ODE:
    #  * ODE: -800\pi u'' + 8\pi u = 0, 0 < x < 100
    #  * BC: u(0) = 10, u(100) = \frac{10}{\cosh(10)}

    [:lin_fin_elt].each do |mthd|
      context "when method: #{mthd}" do
        subject(:numsol) do
          Spitzy::Bvp.new(xrange: [0.0, 100.0], mx: 100, bc: [10.0, 0.00090799859], 
                          a: 800.0*Math::PI, b: 0.0, c: 8.0*Math::PI, f: 0.0, method: mthd)
        end

        it "reports the utilized method" do
          expect(numsol.method).to eq(mthd)
        end

        it "computes a precise numerical solution" do
          exactsol = proc { |x| 10.0 * Math::cosh((100.0 - x)/10.0) / Math::cosh(10.0) }
          u_exact = []
          numsol.x.each { |x| u_exact << exactsol.call(x) }
          error = []
          numsol.u.each_with_index { |pt, i| error << (pt - u_exact[i]).abs }
          maxerror = error.max

          expect(maxerror).to be_within(0.1).of(0.0)
        end
      end 
    end
  end

  describe Spitzy::Ode do

    #  solve the following initial value ODE:
    #  * ODE: y' = -2xy, 0 <= x <= 4,
    #  * IC: y(0) = 1
    #
    #  The exact solution is y = exp(-x^2)

    [:dopri, :euler, :ab2].each do |mthd|
      context "when method: #{mthd}" do
        let(:tol) { 0.1 }

        subject(:numsol) do
          f = proc { |t,y| -2.0*t*y }
          Spitzy::Ode.new(xrange: [0.0,4.0], dx: 0.1, tol: tol,
                          yini: 1.0, method: mthd, &f)
        end

        it "reports the utilized method" do
          expect(numsol.method).to eq(mthd)
        end

        it "computes a precise numerical solution" do
          exact_sol = numsol.x.map { |tt| Math::exp(-(tt**2)) }
          maxerror = exact_sol.each_with_index.map {|n,i| (n - numsol.u[i]).abs }.max

          expect(maxerror).to be_within(tol).of(0.0)
        end
      end 
    end
  end
end
