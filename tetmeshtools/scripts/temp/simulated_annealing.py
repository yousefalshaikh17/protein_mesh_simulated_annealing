import argparse
import sys
import mrcfile
import pathlib
import tetmeshtools.mrclininterpolate as mi
from tetmeshtools.scripts.temp.copy_mesh import read_tetgen_file, remove_suffix, write_tetgen_file
from tetmeshtools.scripts.temp.compute_densities import compute_node_densities, get_density_at_point
from tetmeshtools.meshtools.tetgenstructs import NodePoint
from tetmeshtools.vector3 import Vector3
import tetmeshtools.tetprops as tp
from tetmeshtools.meshtools.trisurface import TriSurface
import numpy as np
import random
from random import Random
import time
import math
import pandas

def get_args():
    description = "Performs simulated annealing given an MRC file and a tetmesh file."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-m",
                        "--mrcfile",
                        type=pathlib.Path,
                        required=True,
                        help="input MRC file")
    
    parser.add_argument("-i",
                        "--tetmeshfile",
                        type=pathlib.Path,
                        required=True,
                        help="input tetmesh file")
    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file name root, (no suffix)")
    parser.add_argument("--density_weight",
                        "-A",
                        type=float,
                        default=1,
                        required=False,
                        help="Sum of Density Weights")
    parser.add_argument("--surface_area_weights",
                        "-B",
                        type=float,
                        default=1,
                        required=False,
                        help="Weights for variance of triangle surface area")
    parser.add_argument("--mutation_probability",
                        type=float,
                        required=True,
                        help="Multiplier controls how large mutation changes can be.")
    parser.add_argument("--mutation_multiplier",
                        type=float,
                        required=True,
                        help="Multiplier controls how large mutation changes can be.")
    parser.add_argument("-t",
                        "--temperature",
                        "--initial_temperature",
                        type=int,
                        required=True,
                        help="Controls the initial temperature.")
    parser.add_argument("--iterations_per_temperature",
                        type=int,
                        required=True,
                        help="How many iterations are done at each temperature change.")
    parser.add_argument("-s",
                        "--seed",
                        type=int,
                        required=False,
                        help="Seed to control randomness.")
    parser.add_argument("--target_density_sum",
                        type=float,
                        required=True,
                        help="Target density sum during the annealing process.")
    parser.add_argument("-k",
                        type=float,
                        required=False,
                        default=0.15,
                        help="Multiplier against T for checking probability of keeping worse fitness.")
    
    
    # Maybe other arguments for mesh optimization later on.
    return parser.parse_args()


def sum_error(target, new):
    return abs(target - new)


def fitness_function(node_densities, surface_area_values, target_density=1, A=1, B=1):
    sum_of_density_error = 0
    for density in node_densities.values():
        sum_of_density_error += abs(target_density - density)

    # print(sum_of_density_error)

    variance = np.var(surface_area_values)
    # print("Raw Density Error and Surface Area Variance: ",sum_of_density_error, variance)
    # print("Weighted Density Error and Surface Area Variance: ",(A * sum_of_density_error),(B * variance))
    # balance = (A * sum_of_density_error)/(B * variance)
    return (A * sum_of_density_error) + (B * variance)#, sum_of_density_error, surface_area_values
        

def mutate(node, temperature, mutation_multiplier=1):
    mutation_max_value = mutation_multiplier * temperature
    new_node = NodePoint(node.index,
        node.x + random.uniform(-mutation_max_value, mutation_max_value),
        node.y + random.uniform(-mutation_max_value, mutation_max_value),
        node.z + random.uniform(-mutation_max_value, mutation_max_value)
    )
    return new_node # Maybe also return how close it was to max value (might help with skipping temperature for faster processing)

def compute_probability_of_keeping(old_fitness, new_fitness, temperature, k=0.15):
    return math.exp((old_fitness - new_fitness)/ (temperature * k))

def simulated_annealing(
        mrc_image,
        nodes,
        faces,
        A=1,
        B=1,
        mutation_probability=0.6,
        target_density=0.6,
        mutation_multiplier = 0.6,
        iterations_per_temp = 100,
        initial_temperature = 100,
        seed = None,
        k = 0.15
    ): # Consider turning into an object

    _early_terminate = 500
    random.seed(seed)
    # Compute initial density values
    node_densities = compute_node_densities(mrc_image, nodes)

    # Compute initial surface area values
    triangle_surface_areas = dict()
    tri_surface = TriSurface(nodes, faces)
    for face_index in tri_surface.get_faces().keys():
        face_nodes = tri_surface.get_triangle_nodes(face_index)
        triangle_surface_areas[face_index] = (tp.area_of_triangle(face_nodes))
    best_fitness = fitness_function(node_densities, list(triangle_surface_areas.values()), target_density, A=A, B=B) #fitness_function(mrc_image, nodes, faces, A=A, B=B, target_density=0.6, cache_densities=True, cache_areas=True)

    current_node = None
    print("Initial Fitness: ", best_fitness)
    # Start
    df = pandas.DataFrame({
        "Temperature": [],
        "Avg Surface Area": [],
        "Surface Area Variance": [],
        "Avg Node Density": [],
        "Node Density Variance": [],
    })
    start_time = time.time()
    try:
        for temperature in range(initial_temperature, 1, -1): # Change decrement method
            temp_start_time = time.time()
            print(f"New Temperature: {temperature}")
            last_best_iteration = 0
            for i in range(iterations_per_temp): # Maybe change it to not be a value
                probabilities = []
                successful_mutations = 0
                old_fitness = best_fitness
                for node in nodes.values():
                    if random.random() <= mutation_probability:
                        current_node = node

                        # mutate node
                        new_node = mutate(node, temperature, mutation_multiplier=mutation_multiplier)

                        # Calculate new density
                        new_density = get_density_at_point(mrc_image, new_node.x, new_node.y, new_node.z)
                        if new_density is None:
                            continue

                        # Save node as replacement for original (Might be swapped back later)
                        nodes[node.index] = new_node
                        
                        

                        old_triangle_areas = dict()
                        for face in faces.values():
                            # Check if a face is associated with a node
                            if new_node.index in {face.vert0, face.vert1, face.vert2}:
                                face_nodes = tri_surface.get_triangle_nodes(face.index)
                                area = tp.area_of_triangle(face_nodes)                                    
                                # Cache old triangle area
                                old_triangle_areas[face.index] = triangle_surface_areas[face.index]
                                # Apply new area
                                triangle_surface_areas[face.index] = area
                                if (np.sign(area) != np.sign(old_triangle_areas[face.index])):
                                    print("inverted.")
                                    exit()

                        # Cache old density
                        old_density = node_densities[node.index]
                        # Apply new density
                        node_densities[node.index] = new_density


                        new_fitness = fitness_function(node_densities, list(triangle_surface_areas.values()), target_density, A=A, B=B)
                        # delta_fitness = best_fitness - new_fitness 
                        keep_probability = compute_probability_of_keeping(best_fitness, new_fitness, temperature, k=k)
                        random_decimal = random.random()
                        if best_fitness > new_fitness or random_decimal <= keep_probability:
                            if best_fitness <= new_fitness and random_decimal <= keep_probability:
                                probabilities.append(keep_probability)
                            if node.index == 1:
                                print(node)
                                print(new_node)
                            current_node = None
                            best_fitness = new_fitness
                            successful_mutations += 1
                            # print(f"New Best Fitness: {best_fitness} at iteration {i} (temperature: {temperature}) (Time Since Start: {round(time.time() - start_time, 2)} seconds)")
                            last_best_iteration = i
                        else:
                            nodes[node.index] = node
                            # revert cached density and triangle area
                            node_densities[node.index] = old_density
                            for face_index, area in old_triangle_areas.items():
                                triangle_surface_areas[face_index] = area
                if len(probabilities) > 0:
                    avg_probability = np.array(probabilities).mean()
                    print(f"Average probability of accepting greater fitness: {avg_probability}")
                if successful_mutations > 0:
                    print(f"New Best Fitness: {best_fitness} (Improvement: {old_fitness-best_fitness}) at iteration {i} (temperature: {temperature}) (successful mutations: {successful_mutations}) (Time Since Start: {round(time.time() - start_time, 2)} seconds)")
                if (i - last_best_iteration) > _early_terminate: # Maybe add some weighted randomness
                    print(f"Early stop for temperature {temperature} at iteration {i}")
                    break
            # Log
            area_values = np.array(list(triangle_surface_areas.values()))
            density_values = np.array(list(node_densities.values()))
            df.loc[len(df.index)] = [
                temperature,
                area_values.mean(),
                np.var(area_values),
                density_values.mean(),
                np.var(density_values)
            ]
    except KeyboardInterrupt:
        if current_node is not None:
            # Swap back for safety measures
            nodes[current_node.index] = current_node
        print(f"Early shutdown.")
    print(f"Finished in {round(time.time() - start_time, 2)} seconds.")
    return nodes, df

def command_validation(args):
    # Verify that the inputs exist
    if not args.mrcfile.exists():
        return f"Error: file {args.mrcfile} does not exist!"
    if not args.tetmeshfile.exists():
        return f"Error: file {args.tetmeshfile} does not exist!"
    if args.mutation_probability < 0:
        return f"Error: Mutation probability must be greater than 0."
    if args.mutation_probability > 1:
        return f"Error: Mutation probability must be less than 1."
    if args.iterations_per_temperature < 1:
        return f"Error: Iterations per temperature must be greater than 0."
    if args.temperature < 1:
        return f"Error: Initial temperature must be greater than 0."
    return None

def main():
    # Verify arguments
    args = get_args()
    error_message = command_validation(args)
    if error_message is not None:
        print(error_message, file=sys.stderr)
        return
    
    # load tetgen mesh 
    nodes, faces, tets, _ = read_tetgen_file(pathlib.Path(remove_suffix(str(args.tetmeshfile))))
    path = pathlib.Path(str(args.mrcfile))

    with mrcfile.mmap(path, mode='r+') as mrc:
        mrc_image = mi.MRCImage(mrc)

    # Assign parameters
    A = args.density_weight
    B = args.surface_area_weights
    mutation_probability= args.mutation_probability
    mutation_multiplier = args.mutation_multiplier
    target_density=args.target_density_sum
    iterations_per_temp = args.iterations_per_temperature
    initial_temperature = args.temperature
    seed = args.seed
    k = args.k
    
    # Perform mesh optimization
    nodes, df = simulated_annealing(
        mrc_image,
        nodes,
        faces,
        A=A,
        B=B,
        target_density=target_density,
        mutation_multiplier=mutation_multiplier,
        mutation_probability=mutation_probability,
        initial_temperature=initial_temperature,
        iterations_per_temp= iterations_per_temp,
        seed=seed,
        k=k
    )
    # Save evolution file
    df.to_csv(str(args.output)+'_evolution.csv', index=False)

    write_tetgen_file(pathlib.Path(str(args.output)), nodes, faces, tets)
    print("Saved Model.")

if __name__ == "__main__":
    main()