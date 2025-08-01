# Ant Colony Optimization for Capacitated Electric Vehicle Routing Problem (ACO-CEVRP)

An implementation of the Ant Colony Optimization (ACO) algorithm to solve the Capacitated Electric Vehicle Routing Problem (CEVRP), addressing the challenge of optimizing electric vehicle routes while considering battery capacity constraints.

## Problem Description

The Capacitated Electric Vehicle Routing Problem (CEVRP) is an extension of the classical Vehicle Routing Problem (VRP) that incorporates the unique constraints of electric vehicles:

- **Limited battery capacity** requiring recharging stations
- **Vehicle capacity constraints** for cargo/passengers
- **Energy consumption** based on distance and load

## Algorithm Overview

This project implements the **Ant Colony Optimization (ACO)** metaheuristic algorithm, which is inspired by the foraging behavior of ants. The algorithm uses:

- **Pheromone trails** to guide solution construction
- **Heuristic information** based on distance and energy consumption
- **Local search** for solution improvement
- **Elitist strategy** for pheromone updates

