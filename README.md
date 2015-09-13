# Spitzy

[![Build Status](https://travis-ci.org/agisga/mixed_models.svg?branch=master)](https://travis-ci.org/agisga/mixed_models)

![Spitzy](spitzy.jpg?raw=true "Optional Title")

Spitzy is this cute pomeranian.
Spitzy reads backwards as *yztips*, which translates into:

***Y*our *Z*appy-*T*appy *I*nitial and boundary value *P*artial (and ordinary) differential equation *S*olver**

Now, spitzy is also a small collection of numerical methods for differential equations, written in Ruby.

## Installation

### Stable

Add this line to your application's Gemfile:

```ruby
gem 'spitzy'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install spitzy

### Development

`spitzy` can be installed using the command line (or something similar):

```
git clone https://github.com/agisga/spitzy.git
cd spitzy/
bundle install
bundle exec rake install
```

The automatic tests can be executed with `bundle exec rspec spec`.

## Documentation, Tutorials and Usage Examples

* Technical Report: <http://agisga.github.io/spitzy/documentation/>

* Slide show: <http://agisga.github.io/presentations/spitzy.html#/>

* Blog posts about some of the implemented methods:

  - [Ordinary differential equations](http://agisga.github.io/ODE/)

  - [Two point boundary value problem](http://agisga.github.io/BVP/)

  - [Advection equation](http://agisga.github.io/Advection-Equation/)

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request
