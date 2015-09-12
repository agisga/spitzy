# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'spitzy/version'

Gem::Specification.new do |spec|
  spec.name          = "spitzy"
  spec.version       = Spitzy::VERSION
  spec.authors       = ["Alexej Gossmann"]
  spec.email         = ["alexej.go@googlemail.com"]

  spec.summary       = %q{Solve differential equations in pure Ruby.}
  spec.description   = %q{A toolbox of numerical differential equation solvers written in pure Ruby.}
  spec.homepage      = "https://github.com/agisga/spitzy.git"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.9"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.2"

  spec.add_runtime_dependency "nmatrix", "~> 0.2.0"
  spec.add_runtime_dependency "nmatrix-lapacke", "~> 0.2.0"

  # This gem will work with Ruby 2.1 or newer, because it uses required 
  # keyword arguments (i.e. keyword arguments without default values)
  spec.required_ruby_version = '>= 2.1'
end
