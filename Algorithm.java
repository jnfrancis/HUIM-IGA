import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Comparator;
import java.util.Arrays;
import java.util.function.Predicate;
import java.util.function.Function;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

public class Algorithm {
	double maxMemory = 0; // the maximum memory usage
	long startTimestamp = 0; // the time the algorithm started
	long endTimestamp = 0; // the time the algorithm terminated	
	final int pop_size = 20;//the size of population
	final int max_evalutions = 60000; //the number of function evalutions
	int ev_num = 0; // (Current) evaluation times
	int chrom_len; //Length of chromosome

	Map<Integer, Integer> mapItemToTWU;
	Map<Integer, Integer> mapItemToTWU0;
	List<Integer> twuPattern;// the items which has twu value more than minUtil	
	BufferedWriter writer = null; // writer to write the output file

	
	// this class represent an item and its utility in a transaction
	class Pair {
		int item = 0;
		int utility = 0;
	}
	
	class A_Chrom  implements Cloneable{
		BitSet chromosome;
		int fitness;

		A_Chrom(int chrom_len){
			chromosome = new BitSet(chrom_len);
		}

		@Override
		public A_Chrom clone(){
			A_Chrom a_chrom = null;
			try {
				a_chrom = (A_Chrom)super.clone();
				a_chrom.chromosome = (BitSet) chromosome.clone();
			}catch(CloneNotSupportedException e){
				e.printStackTrace();
			}
			return a_chrom;
		}

	}
	
	class HUI {
		String itemset;
		int fitness;

		public HUI(String itemset, int fitness) {
			super();
			this.itemset = itemset;
			this.fitness = fitness;
		}

	}
	
	//use Item to create bitmap
	class Item
	{
		int item;
		BitSet TIDS;

		public Item(){
			TIDS =  new BitSet(database.size());
		}

		public Item(int item){
			TIDS=new BitSet(database.size());
			this.item=item;	
		}
	}
	
	List<A_Chrom> HUIsChrom = new ArrayList<A_Chrom>();
	List<HUI> huiSets = new ArrayList<HUI>();

	List<List<Pair>> database = new ArrayList<List<Pair>>();
	List<Item> Items; //bitmap database representation
	List<Integer> twu_value = new ArrayList<Integer>();


	/**
	 * Calculating the fitness value of a chromosome
	 */
	public void getFitnessValue(A_Chrom a_chrom, List<Integer> templist)
	{
		BitSet a_chrom_bitset = a_chrom.chromosome;
		int k = a_chrom_bitset.cardinality();
		if (k == 0)
		{
			a_chrom.fitness=0;
			return;
		}

		int i, p, q, temp,m;
		int sum, fitness = 0;
		for (m = 0; m < templist.size(); m++)
		{
			p = templist.get(m).intValue();
			i = 0;
			q = 0;
			temp = 0;
			sum = 0;

			while (q < database.get(p).size() && i < chrom_len)
			{
				if(a_chrom.chromosome.get(i))
				{
					if (database.get(p).get(q).item == twuPattern.get(i))
					{
						sum = sum + database.get(p).get(q).utility;
						++i;
						++q;
						++temp;
					}
					else
					{
						++q;
					}
				}
				else
				{
					++i;
				}
			}
			if (temp == k){
				fitness = fitness + sum;
			}
		}
		a_chrom.fitness = fitness;

		//System.out.println(ev_num+" "+ HUIsChrom.size());
		ev_num++;
	}


	public void runAlgorithm(String input, String output, int minUtility)
			throws IOException {
		// reset maximum
		maxMemory = 0;
		startTimestamp = System.currentTimeMillis();
		writer = new BufferedWriter(new FileWriter(output));
		// We create a map to store the TWU of each item
		mapItemToTWU = new HashMap<Integer, Integer>();
		mapItemToTWU0 = new HashMap<Integer, Integer>();

		// We scan the database a first time to calculate the TWU of each item.
		BufferedReader myInput = null;
		String thisLine;
		try {
			// prepare the object for reading the file
			myInput = new BufferedReader(new InputStreamReader(
					new FileInputStream(new File(input))));
			// for each line (transaction) until the end of file
			while ((thisLine = myInput.readLine()) != null) {
				// if the line is a comment, is empty or is a
				// kind of metadata
				if (thisLine.isEmpty() == true || thisLine.charAt(0)== '#'
						|| thisLine.charAt(0) == '%'
						|| thisLine.charAt(0) == '@') {
					continue;
				}

				// split the transaction according to the : separator
				String split[] = thisLine.split(":");
				// the first part is the list of items
				String items[] = split[0].split(" ");
				// the second part is the transaction utility
				int transactionUtility = Integer.parseInt(split[1]);
				// for each item, we add the transaction utility to its TWU
				for (int i = 0; i < items.length; i++) {
					// convert item to integer
					Integer item = Integer.parseInt(items[i]);
					// get the current TWU of that item
					Integer twu = mapItemToTWU.get(item);
					Integer twu0 = mapItemToTWU0.get(item);
					// add the utility of the item in the current transaction to
					// its twu
					twu = (twu == null) ? transactionUtility : twu
							+ transactionUtility;
					twu0 = (twu0 == null) ? transactionUtility : twu0
							+ transactionUtility;
					mapItemToTWU.put(item, twu);
					mapItemToTWU0.put(item, twu0);
				}
			}
		} catch (Exception e) {
			// catches exception if error while reading the input file
			e.printStackTrace();
		} finally {
			if (myInput != null) {
				myInput.close();
			}
		}
		// SECOND DATABASE PASS TO CONSTRUCT THE DATABASE
		// OF 1-ITEMSETS HAVING TWU >= minutil (promising items)
		try {
			// prepare object for reading the file
			myInput = new BufferedReader(new InputStreamReader(
					new FileInputStream(new File(input))));
			// variable to count the number of transaction
			// for each line (transaction) until the end of file
			while ((thisLine = myInput.readLine()) != null) {
				// if the line is a comment, is empty or is a
				// kind of metadata
				if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#'
						|| thisLine.charAt(0) == '%'
						|| thisLine.charAt(0) == '@') {
					continue;
				}

				// split the line according to the separator
				String split[] = thisLine.split(":");
				// get the list of items
				String items[] = split[0].split(" ");
				// get the list of utility values corresponding to each item
				// for that transaction
				String utilityValues[] = split[2].split(" ");

				// Create a list to store items and its utility
				List<Pair> revisedTransaction = new ArrayList<Pair>();
				// Create a list to store items
				List<Integer> pattern = new ArrayList<Integer>();
				// for each item
				for (int i = 0; i < items.length; i++) {
					// / convert values to integers
					Pair pair = new Pair();
					pair.item = Integer.parseInt(items[i]);
					pair.utility = Integer.parseInt(utilityValues[i]);
					// if the item has enough utility
					if (mapItemToTWU.get(pair.item) >= minUtility) {
						// add it
						revisedTransaction.add(pair);
						pattern.add(pair.item);
					}else{
						mapItemToTWU0.remove(pair.item);
					}
				}
				// Copy the transaction into database but
				// without items with TWU < minutility
				database.add(revisedTransaction);
			}
		} catch (Exception e) {
			// to catch error while reading the input file
			e.printStackTrace();
		} finally {
			if (myInput != null) {
				myInput.close();
			}
		}
		
		twuPattern = new ArrayList<Integer>(mapItemToTWU0.keySet());
		Collections.sort(twuPattern); //
		Items = new ArrayList<Item>();
		for(Integer tempitem:twuPattern){
			Items.add(new Item(tempitem.intValue()));
		}
	
		//scan database to create bitmap
		for(int i=0;i<database.size();++i){
			for(int j=0;j<Items.size();++j){
				for(int k=0;k<database.get(i).size();++k){
					if(Items.get(j).item==database.get(i).get(k).item){
						Items.get(j).TIDS.set(i);
					}
				}
			}
		}
		// check the memory usage
		checkMemory();

		// Mine the database recursively
		if (twuPattern.size() > 0) {
			chrom_len = twuPattern.size();
			int tempA,tempB,tempA1,tempB1;
			int temp1 = 0, temp2 = 0;
			int choosed_HUI_index_1, choosed_HUI_index_2;
			int HUI_pointer = 0;

			for (Integer key:mapItemToTWU0.keySet())
			{
				twu_value.add(mapItemToTWU0.get(key));
			}

			Integer[] TWU_index = new Integer[chrom_len];
			int ind;
			for(ind=0;ind<chrom_len;ind++)
				TWU_index[ind]=ind;

			Arrays.sort(TWU_index, new Comparator<Integer>() {
				@Override public int compare(final Integer i, final Integer j) {
					return -Float.compare(twu_value.get(i), twu_value.get(j)); //Descending sort
				}
			});

			List<A_Chrom> population = new ArrayList<A_Chrom>();
			List<A_Chrom> subPopulation = new ArrayList<A_Chrom>();

			// initial population
			Initialization(minUtility,TWU_index,population);

			int m,j;
			int sum_fitness_value = 0;
			int[] fitness_values = new int[pop_size];

			while (ev_num < max_evalutions)
			{
				if(HUIsChrom.size() == 1)
				{
					tempA1 = (int)(Math.random()*pop_size);
					tempB1 = (int)(Math.random()*pop_size);

					A_Chrom replace_indiv1 = HUIsChrom.get(0).clone();
					population.set(tempA1,replace_indiv1);
					population.set(tempB1,replace_indiv1);
				}

				if(HUIsChrom.size() > 1)
				{
					tempA = (int)(Math.random()*pop_size);
					tempB = (int)(Math.random()*pop_size);
					if (HUIsChrom.size() == 2)
					{
						choosed_HUI_index_1 = 0;
						choosed_HUI_index_2 = 1;
					}
					else
					{
						if (HUI_pointer + 2 > HUIsChrom.size()-1)
						{
							HUI_pointer = 0;
							choosed_HUI_index_1 = HUI_pointer;
							choosed_HUI_index_2 = HUI_pointer + 1;
						}
						else
						{
							HUI_pointer = HUI_pointer + 2;
							if (HUI_pointer + 1 > HUIsChrom.size()-1)
							{
								choosed_HUI_index_1 = HUI_pointer;
								choosed_HUI_index_2 = 0;
								HUI_pointer = -1;
							}
							else {
								choosed_HUI_index_1 = HUI_pointer;
								choosed_HUI_index_2 = HUI_pointer + 1;
							}
						}
					}

					A_Chrom replace_choose_HUI_1 = HUIsChrom.get(choosed_HUI_index_1).clone();
					A_Chrom replace_choose_HUI_2 = HUIsChrom.get(choosed_HUI_index_2).clone();
					// replace
					population.set(tempA,replace_choose_HUI_1);
					population.set(tempB,replace_choose_HUI_2);
				}

				List<Integer> accumulative_fitness_value = new ArrayList<Integer>();

				m=0;
				for(A_Chrom ac1:population)
					fitness_values[m++] = ac1.fitness;
				sum_fitness_value  = totalSum_and_accumulative(fitness_values, accumulative_fitness_value);

				while (subPopulation.size() < pop_size)
				{
					// roulette wheel selection
					temp1 = rouletteSelection(sum_fitness_value,accumulative_fitness_value);
					temp2 = rouletteSelection(sum_fitness_value,accumulative_fitness_value);

					while (temp1 == temp2) {
						temp2=(temp2+(int)(Math.random()*1000))%pop_size;
					}
					uniform_crossover(population, subPopulation, temp1, temp2, TWU_index,minUtility);
				}
				// mutation
				Mutation(subPopulation, TWU_index, minUtility);

				// Merged population
				subPopulation.addAll(population);

				// Delete duplicate individuals
				List<A_Chrom> distinct_pop = subPopulation.stream()
						.filter( distinctByChromosome(p -> p.chromosome) )
						.collect( Collectors.toList() );

				Collections.sort(distinct_pop,new Comparator<A_Chrom>(){
					public int compare(A_Chrom arg0, A_Chrom arg1) {
						return arg1.fitness - arg0.fitness ;
					}
				});

                if (distinct_pop.size() < pop_size)
                {
                    int need_num = pop_size - distinct_pop.size();
                    for (int nedd_i = 0; nedd_i< need_num; nedd_i++)
                        distinct_pop.add(distinct_pop.get(nedd_i));
                }

				for (j = 0; j < pop_size; j++) {
					population.set(j, distinct_pop.get(j));
				}
				subPopulation.clear();
			}
		}

        HUIChroms_to_HUIs(); // Decode

		writeOut();
		// check the memory usage again and close the file.
		checkMemory();
		// close output file
		writer.close();
		// record end time
		endTimestamp = System.currentTimeMillis();

	}



	private void Initialization(int minUtility, Integer[] TWU_index, List<A_Chrom> population)
	{
		int i = 0, j, k, choosed_item;
		List<Integer> transList;
		List<Integer> accumulative_twu_value = new ArrayList<Integer>();
		//Integer sum_twu_value  = totalSum_and_accumulative(twu_value, accumulative_twu_value);

		int[] twu_arr = twu_value.stream().mapToInt(Integer::valueOf).toArray();
		int sum_twu_value = totalSum_and_accumulative(twu_arr, accumulative_twu_value);


		while (i < pop_size)
		{
			A_Chrom tempChrom = new A_Chrom(chrom_len);
			j = 0;
			k = (int) (Math.random() * chrom_len);

			while (k == 0)
			{
				k = (int) (Math.random() * chrom_len);
			}

			while (j < k)
			{
				choosed_item = rouletteSelection(sum_twu_value,accumulative_twu_value);

				if (!tempChrom.chromosome.get(choosed_item)) {
					j++;
					tempChrom.chromosome.set(choosed_item);
				}
			}

			transList=new ArrayList<Integer>();
			Indiv_Repair(tempChrom, transList, TWU_index);
			getFitnessValue(tempChrom, transList);
			population.add(tempChrom);

			if (tempChrom.fitness >= minUtility&&tempChrom.chromosome.cardinality()>0) {
				addHUI(tempChrom);
			}
			i++;
		}
	}


	public static <T> Predicate<T> distinctByChromosome(Function<? super T, Object> keyExtractor)
	{
		Map<Object, Boolean> map = new ConcurrentHashMap<>();
		return t -> map.putIfAbsent(keyExtractor.apply(t), Boolean.TRUE) == null;
	}


	private int totalSum_and_accumulative(int[] value, List<Integer> accumulative_value)
	{
		int tempSum = 0;
		for (int a:value)
		{
			tempSum += a;
			accumulative_value.add(tempSum);
		}
		return tempSum;
	}


	private int rouletteSelection(int sum, List<Integer> accumulative_value)
	{
		double roulette_position = sum * Math.random();
		return BinarySearch(accumulative_value, roulette_position);
	}

	private int BinarySearch(List<Integer> accumulative_value, double key){
		int low = 0, high = accumulative_value.size()-1, mid;
		if (key > accumulative_value.get(high))
			return high;
		if (key < accumulative_value.get(low))
			return low;
		while (low <= high){
			mid = (low + high) / 2;
			if (accumulative_value.get(mid) == key)
				return mid;
			if (accumulative_value.get(mid) < key && accumulative_value.get(mid+1) >= key)
				return mid + 1;
			if (accumulative_value.get(mid) > key)
				high = mid - 1;
			if (accumulative_value.get(mid) < key)
				low = mid + 1;
		}
		System.out.println("something wrong!");
		return -1;
	}


	/**
	 * Uniform crossover operator
	 */
	private void uniform_crossover(List<A_Chrom> population, List<A_Chrom> subPopulation, int temp1, int temp2, Integer[] TWU_index, int minUtility) {
		int j,i;
		double r;
		boolean temp_p1,temp_p2;
		List<Integer> transList_1;
		List<Integer> transList_2;

		A_Chrom temp1Chrom = population.get(temp1).clone();
		A_Chrom temp2Chrom = population.get(temp2).clone();

		for(j = 0; j < chrom_len; j++)
		{
			r = Math.random();
			if (r>0.5)
			{
				temp_p1 = temp1Chrom.chromosome.get(j);
				temp_p2 = temp2Chrom.chromosome.get(j);
				if (temp_p1 == temp_p2)
					continue;
				else
				{
					if (temp_p1)
					{
						temp1Chrom.chromosome.clear(j);
						temp2Chrom.chromosome.set(j);
					}
					else
					{
						temp2Chrom.chromosome.clear(j);
						temp1Chrom.chromosome.set(j);
					}
				}
			}
		}

		transList_1=new ArrayList<Integer>();
		Indiv_Repair(temp1Chrom, transList_1, TWU_index);

		if(HUIsChrom.size() != 0)
		{
			for(i = 0; i < HUIsChrom.size(); ++i)
			{
				if(temp1Chrom.chromosome.equals(HUIsChrom.get(i).chromosome))
				{
					subPopulation.add(temp1Chrom);
					break;
				}
				if (i == (HUIsChrom.size() -1))
				{
					getFitnessValue(temp1Chrom, transList_1);
					subPopulation.add(temp1Chrom);
					if (temp1Chrom.fitness >= minUtility&&temp1Chrom.chromosome.cardinality()>0) {
						A_Chrom temp_pa1 = temp1Chrom.clone();
						HUIsChrom.add(temp_pa1);
					}
				}
			}
		}
		else
		{
			getFitnessValue(temp1Chrom, transList_1);
			subPopulation.add(temp1Chrom);
			if (temp1Chrom.fitness >= minUtility&&temp1Chrom.chromosome.cardinality()>0) {
				A_Chrom temp_pa1 = temp1Chrom.clone();
				HUIsChrom.add(temp_pa1);
			}
		}


		transList_2=new ArrayList<Integer>();
		Indiv_Repair(temp2Chrom, transList_2, TWU_index);

		if(HUIsChrom.size() != 0)
		{
			for(i = 0; i < HUIsChrom.size(); ++i){

				if(temp2Chrom.chromosome.equals(HUIsChrom.get(i).chromosome)){
					subPopulation.add(temp2Chrom);
					break;
				}

				if (i == (HUIsChrom.size() -1))
				{
					getFitnessValue(temp2Chrom, transList_2);
					subPopulation.add(temp2Chrom);
					if (temp2Chrom.fitness >= minUtility&&temp2Chrom.chromosome.cardinality()>0) {
						A_Chrom temp_pa2 = temp2Chrom.clone();
						HUIsChrom.add(temp_pa2);
					}
				}
			}
		}
		else
		{
			getFitnessValue(temp2Chrom, transList_2);
			subPopulation.add(temp2Chrom);
			if (temp2Chrom.fitness >= minUtility&&temp2Chrom.chromosome.cardinality()>0) {
				A_Chrom temp_pa2 = temp2Chrom.clone();
				HUIsChrom.add(temp_pa2);
			}
		}

	}


	/**
	 * Mutation operator
	 */
	private void Mutation(List<A_Chrom> subPopulation, Integer[] TWU_index, int minUtility) {
        int i;
		List<Integer> transList;

		for (i = 0; i < pop_size; i++)
		{
			// single point mutation (one gene is always randomly selected and flipped)
			int temp = (int) (Math.random() * chrom_len);
			if (subPopulation.get(i).chromosome.get(temp)) {
				subPopulation.get(i).chromosome.clear(temp);
			} else {
				subPopulation.get(i).chromosome.set(temp);
			}
			transList=new ArrayList<Integer>();
			Indiv_Repair(subPopulation.get(i), transList, TWU_index);

			if(HUIsChrom.size() != 0)
			{
			    int k;
				for(k = 0; k < HUIsChrom.size(); k++)
				{
					if(subPopulation.get(i).chromosome.equals(HUIsChrom.get(k).chromosome))
					{
						HUIsNE(subPopulation.get(i));
						transList=new ArrayList<Integer>();
						Indiv_Repair(subPopulation.get(i), transList, TWU_index);
						break;
					}
				}
			}
			getFitnessValue(subPopulation.get(i), transList);
			// insert chromosome has higher utility into huiSets
			if (subPopulation.get(i).fitness >= minUtility&&subPopulation.get(i).chromosome.cardinality()>0) {
				addHUI(subPopulation.get(i));
			}
		}
	}



	/**
	 * neighborhood exploration
	 */
	public void HUIsNE(A_Chrom now_indiv)
	{
		BitSet temp_nowindiv = (BitSet)now_indiv.chromosome.clone();
		int point1 = (int)(chrom_len*Math.random());
		if (now_indiv.chromosome.get(point1))
		{
			now_indiv.chromosome.clear(point1);
			now_indiv.chromosome.set(temp_nowindiv.nextClearBit(0));
		}
		else
		{
			now_indiv.chromosome.set(point1);
			now_indiv.chromosome.clear(temp_nowindiv.nextSetBit(0));
		}
	}


	/**
	 * Individual repair
	 */
	public void Indiv_Repair(A_Chrom now_indiv, List<Integer> transList, Integer[] TWU_index)
	{
		int i,j;

		if(now_indiv.chromosome.isEmpty()){
			return ;
		}
		BitSet tempBitSet = new BitSet(database.size());
		BitSet midBitSet = (BitSet)tempBitSet.clone();

		for (i = 0; i < chrom_len; i++)
		{
			if (now_indiv.chromosome.get(TWU_index[i]))
			{
				tempBitSet = (BitSet)Items.get(TWU_index[i]).TIDS.clone();
				break;
			}
		}

		for (j = i + 1; j < chrom_len; j++)
		{
			if (now_indiv.chromosome.get(TWU_index[j]))
			{
				tempBitSet.and(Items.get(TWU_index[j]).TIDS);

				if(tempBitSet.cardinality() != 0)
				{
					midBitSet = (BitSet)tempBitSet.clone();
				}
				else{
					tempBitSet = (BitSet)midBitSet.clone();
					now_indiv.chromosome.clear(TWU_index[j]);
				}
			}
		}

		if (tempBitSet.cardinality() != 0)
		{
			int m;
			for (m=0;m<tempBitSet.length();m++)
			{
				if (tempBitSet.get(m))
					transList.add(m);
			}
		}
	}



	/**
	 * add hui
	 */
	private void addHUI(A_Chrom tempBAIndividual){
		int i;
	    if(HUIsChrom.size() != 0)
			for (i=0; i< HUIsChrom.size(); i++)
				if(tempBAIndividual.chromosome.equals(HUIsChrom.get(i).chromosome))
					return;

		A_Chrom temp_indiv = tempBAIndividual.clone();
		HUIsChrom.add(temp_indiv);
	}


    /**
     * Decode
     */
    private void HUIChroms_to_HUIs()
    {
        int k,i;
        A_Chrom temp_chromosome;
        for (k = 0; k < HUIsChrom.size(); k++)
        {
            StringBuilder temp = new StringBuilder();
            temp_chromosome = HUIsChrom.get(k);
            for (i = 0; i < twuPattern.size(); i++)
            {
                if (temp_chromosome.chromosome.get(i)) {
                    temp.append(twuPattern.get(i));
                    temp.append(' ');
                }
            }
            huiSets.add(new HUI(temp.toString(), temp_chromosome.fitness));
        }
    }


	/**
	 * Method to write a high utility itemset to the output file.
	 * 
	 * @throws IOException
	 */
	private void writeOut() throws IOException {
		// Create a string buffer
		StringBuilder buffer = new StringBuilder();

//        for (int i = 0; i < save_res.size(); i++) {
//            buffer.append(i);
//            buffer.append(" ");
//            buffer.append(save_res.get(i).get(i+1));
//            if(i != save_res.size() -1){
//                buffer.append(System.lineSeparator());
//            }
//        }


		// append the prefix
		for (int i = 0; i < huiSets.size(); i++) {
			buffer.append(huiSets.get(i).itemset);
			// append the utility value
			buffer.append("#UTIL: ");
			buffer.append(huiSets.get(i).fitness);
			if(i != huiSets.size() -1){
				buffer.append(System.lineSeparator());
			}
		}

		// write to file
		writer.write(buffer.toString());
	}
	
	/**
	 * Method to check the memory usage and keep the maximum memory usage.
	 */
	private void checkMemory() {
		// get the current memory usage
		double currentMemory = (Runtime.getRuntime().totalMemory() - Runtime
				.getRuntime().freeMemory()) / 1024d / 1024d;
		// if higher than the maximum until now
		if (currentMemory > maxMemory) {
			// replace the maximum with the current memory usage
			maxMemory = currentMemory;
		}
	}
	
	/**
	 * Print statistics about the latest execution to System.out.
	 */
	public void printStats() {
		System.out
				.println("=============  HUIM-IGA ALGORITHM  =============");
		System.out.println(" Total time ~ " + (endTimestamp - startTimestamp)
				+ " ms");
		System.out.println(" Memory ~ " + maxMemory + " MB");
		System.out.println(" High-utility itemsets count : " + huiSets.size());
		System.out
				.println("===================================================");
	}
	

}
