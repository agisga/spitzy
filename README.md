# Spitzy

![Spitzy](spitzy.jpg?raw=true "Optional Title")

Spitzy is this cute pomeranian.
Spitzy reads backwards as *yztips*, which translates into:

***Y*our *Z*appy-*T*appy *I*nitial and boundary value *P*artial (and ordinary) differential equation *S*olver**

Now, spitzy is also a growing collection of numerical methods for differential equations, written in Ruby.
To my knowledge, apart from an [interface with the DASPK Fortran library](https://rubygems.org/gems/rb-daspk/versions/0.0.7-x86-mswin32-60), there currently does not exist another differential equation solver gem for Ruby.

<!--
Welcome to your new gem! In this directory, you'll find the files you need to be able to package up your Ruby library into a gem. Put your Ruby code in the file `lib/spitzy`. To experiment with that code, run `bin/console` for an interactive prompt.

TODO: Delete this and the text above, and describe your gem 
-->

## Installation

<!--
Add this line to your application's Gemfile:

```ruby
gem 'spitzy'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install spitzy
-->

Ruby is required in version >=2.0 because keyword arguments are excessively used in `spitzy`.
Moreover, prior to the installation of `spitzy`, currently the `NMatrix` gem needs to be installed in its development version (because `Poissons_eq` uses `#meshgrid`)
 from <https://github.com/SciRuby/nmatrix.git>.

Then `spitzy` can be installed using the command line (or something similar):

```
git clone https://github.com/agisga/spitzy.git
cd spitzy/
rake install
```

Feel free to contact me at alexej.go [at] googlemail.com, in case of difficulties with the installation.

## Documentation and Usage Examples

* Documentation: <http://agisga.github.io/spitzy/documentation/>

* Slide show: <http://agisga.github.io/presentations/spitzy.html#/>

* Blog posts about some of the implemented methods:

  - [Ordinary differential equations](http://agisga.github.io/ODE/)

  - [Two point boundary value problem](http://agisga.github.io/BVP/)

  - [Advection equation](http://agisga.github.io/Advection-Equation/)


<!--
## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release` to create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).
-->

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request
