---
layout: post
title: Monte Carlo simulation of airline overbooking
tags: [optimization]
snippet: Airlines overbook flights to maximize profits. But by how much should they overbook? Simulations can address this question!
author: Mira Khare, Melanie Huynh, Arni Sturluson, Cory Simon
---

Airlines [overbook](https://en.wikipedia.org/wiki/Overselling) flights to maximize profits.

Consider a flight with seats for 100 passengers.

If the airline allows only up to 100 reservations for the flight, each customer will be guarenteed a seat on the plane. However, several seats will likely be empty because a fraction of the customers that reserve a seat will likely miss the flight. For each empty seat, the airline looses revenue.

If the airline allows greater than 100 reservations-- i.e. if the airline overbooks the flight-- it is more likely that enough customers show up to fill the plane so that the airline receives revenue from each seat on the airplane. However, if more than 100 customers that reserved a flight show up, then (i) customers bumped from the flight will be angry and not fly with that airline again and (ii) the airline must pay out vouchers to incentivize volunteers to take a different flight. This is costly.

Clearly, there is a delicate balance of how much to overbook. Don't overbook at all: all customers are happy, but flights operate below full capacity and we lose out on revenue. Overbook too much: flights are likely operating near full capacity so we receive revenue from each seat, but customers bumped from the flight are angry and costly vouchers must be paid to volunteers.

We wrote a code in the [Julia programming language](https://julialang.org/) to simulate the process of customers showing up to flights and identify the optimal amount to overbook.

As input, we need the probability that a given customer shows up to the flight after making a reservation.

```julia
probability_show = 0.935
```

The function `show_up` simulates the stochastic process of a customer showing up for the flight after making a reservation.

```julia
function show_up(probability_show::Float64)
  if rand() <= probability_show
    return true # passenger showed up
  else
    return false # passenger didn't show up
  end
end
```

The function `simulate_flight` simulates the stochastic process of `nb_tickets_sold` reservations made for a flight and returns the number of customers that show up.

```julia
function simulate_flight(nb_tickets_sold::Int64, 
                         probability_show::Float64)
  n = 0 # number of folks who bought tix that will show up
  for i = 1:nb_tickets_sold
    if show_up(probability_show)
      n = n + 1
    end
  end
  return n
end
```

The revenue received from a flight will depend on the number of reservations made (`nb_tickets_sold`), number of seats on the plane (`nb_seats`), the probability that a given customer shows up (`probability_show`), the revenue received per seat (`revenue_per_seat`), and the cost of a voucher to incentivize volunteers to take a different flight when the plane is overbooked, (`voucher_cost`). The latter may include the implicit cost of losing customers to different airlines after angering them. The function `simulate_net_revenue`simulates the stochastic process and returns the net revenue from the flight.

```julia
function simulate_net_revenue(nb_tickets_sold::Int64, 
                              nb_seats::Int, 
                              probability_show::Float64, 
                              revenue_per_seat::Float64, 
                              voucher_cost::Float64)
  # how many ticket purchasers actually showed up?
  nb_shows = simulate_flight(nb_tickets_sold, probability_show)
  # no one bumped from flight if less or equal folks show up than for
  #   the number of seats we have
  if nb_shows <= nb_seats
    return revenue_per_seat * nb_shows
  # if more customers show up than seats we hv, must pay out vouchers
  else
    vouchers_out = nb_shows - nb_seats
    return nb_seats * revenue_per_seat - voucher_cost * vouchers_out
  end
end
```

For a $350 flight on a plane of 100 seats, where vouchers are twice the cost of the ticket:

```julia
nb_seats = 100

revenue_per_seat = 350.0 # USD
voucher_cost = revenue_per_seat * 2.0 # USD
```

We now simulate the process of different amounts of overbooking to identify the optimal amount of overbooking. Since this is a stochastic process, we run 10,000 simulations to simulate 10,000 flights so we can calculate the average net revenue for each amount of overbooking and gauge the variance among flights.

```julia
nb_flights = 10000
max_overbooking = 15

# revenue[i, k] is net revenue received from flight i with k-1 tickets over capacity sold
revenue = zeros(nb_flights, max_overbooking + 1)
for tix_overbooked = 0:max_overbooking
  nb_tickets_sold = nb_seats + tix_overbooked
  for f = 1:nb_flights # simulate nb_flights flights
    revenue[f, tix_overbooked + 1] =
        simulate_net_revenue(nb_tickets_sold, 
                             nb_seats, 
                             probability_show, 
                             revenue_per_seat, 
                             voucher_cost)
  end
end
```

We depict the results using a box plot for each amount of overbooking:

```julia
using PyPlot

boxplot(revenue, labels=0:max_overbooking)
xlabel("# tickets sold beyond capacity")
ylabel("net revenue")
```

<figure>
    <img src="/images/overbooking_boxplot.png" alt="image" style="width: 100%;">
    <figcaption>Box plot showing the spread of net revenue among 10,000 simulated flights depending on the number of tickets sold beyond capacity.
    </figcaption>
</figure>

The box plot shows that, if we do not overbook (0 on the x-axis), we on average receive less net revenue than if we do overbook. If we overbook too much, e.g. sell 15 tickets beyond capacity, we see that the average net revenue is less than if we overbook because we are paying out costly vouchers. The optimal amount of overbooking with these parameters is shown by the box plot to be 5-7 tickets beyond capacity. Of course, this result depends on the probability that a customer will show up, the cost of the flight and voucher, and number of seats on the plane.

Note: airlines hire analysts to do exactly this. See [this NYT article](https://www.nytimes.com/2007/05/30/business/30bump.html?pagewanted=all&_r=0).

This code was developed by Mira Khare, Melanie Huynh, Arni Sturluson, and two high school students for the [Summer Experience in Science and Engineering for Youth](http://cbee.oregonstate.edu/sesey) (SESEY) at Oregon State University.
