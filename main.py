import random

# get random num in range
def get_random(array, num):
    return random.sample(array, num)

# get random gene
def get_gene(gene_possible_array, gene_length):
    return get_random(gene_possible_array, gene_length) + [-1]

# get a gene array
def get_random_genes(population, gene_possible_array, gene_length):
    gene_array = []
    for i in range(population):
        gene_array.append(get_gene(gene_possible_array, gene_length))
    return gene_array


class ChessboardGene:
    def __init__(self, population: int, gene_possible_array: list, gene_length: int):
        self.gene_array = get_random_genes(population, gene_possible_array, gene_length)
        self.gene_possible_array = gene_possible_array
        self.gene_length = gene_length
        self.min_cost = -1
        self.population = len(self.gene_array)

    # constructor
    def generate_genes(self):
        self.gene_array = get_random_genes(self.population, self.gene_possible_array, self.gene_length)
    
    # new generation
    def cross_over(self, chosen_rate=0.4):
        split_index = random.randint(0, self.gene_length)
        child_genes = []
        for i in range(0, len(self.gene_array[:int(len(self.gene_array) * chosen_rate)]), 2):
            if len(self.gene_array) <= i + 1:
                break
            father_part1 = self.gene_array[i][:split_index]
            father_part2 = self.gene_array[i][split_index:]
            mother_part1 = self.gene_array[i + 1][:split_index]
            mother_part2 = self.gene_array[i + 1][split_index:]
            first_child = father_part1 + mother_part2
            second_child = mother_part1 + father_part2
            child_genes.append(first_child)
            child_genes.append(second_child)
        self.gene_array = self.gene_array + child_genes
    
    # mutation some genes
    def mutation(self):
        all_genes = self.gene_array
        mutated_genes = []
        for gene in all_genes:
            mutated_gene = gene
            rand = random.randint(1, self.gene_length - 1)
            mutated_gene[rand] = get_random(self.gene_possible_array, 1)[0]
            mutated_genes.append(mutated_gene)
        self.gene_array = self.gene_array + mutated_genes

    # evaluate generation
    def evaluation(self, chosen_rate=0.6):
        genes_threats = []
        for gene in self.gene_array:
            threat = 0
            for i in range(0, len(gene) - 1):
                for j in range(0, len(gene) - 1):
                    if i == j:
                        break
                    if gene[i] == gene[j]:
                        threat = threat + 1
                    if abs(gene[i] - gene[j]) == abs(i - j):
                        threat = threat + 1
            gene[8] = threat
            genes_threats.append(gene)
        top_genes = sorted(genes_threats, key=lambda items: items[8])[:int(len(self.gene_array) * chosen_rate)]
        self.min_cost = top_genes[0][8]
        self.gene_array = top_genes

    # run functions in loop
    def run(self, top=3):
        i = 0
        while self.min_cost != 0:
            i = i + 1
            self.cross_over()
            self.mutation()
            self.evaluation()
            print(i, ")", self.gene_array[:top])


genes = ChessboardGene(10, range(1, 9), 8)
genes.run(7)
