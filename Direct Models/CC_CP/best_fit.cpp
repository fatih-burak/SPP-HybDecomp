#include "model.h"
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
using namespace std;

// Finds the lowest "niche" in the skyline
// Returns: (min_height_value, start_index, gap_width)
tuple<int, int, int> findLowestGap(const vector<int>& skyline_array) {
	// 1. Find the minimum height value
	int minVal = *min_element(skyline_array.begin(), skyline_array.end());

	// 2. Find the first index with this min value
	int n = skyline_array.size();
	int startIdx = -1;
	for (int i = 0; i < n; i++) {
		if (skyline_array[i] == minVal) {
			startIdx = i;
			break;
		}
	}
	if (startIdx == -1) return { -1, -1, 0 }; // should never happen

	// 3. Count how many consecutive cells have the same min value
	int width = 0;
	for (int i = startIdx; i < n && skyline_array[i] == minVal; i++) {
		width++;
	}

	//cout << "Lowest height = " << minVal << "\n";
	//cout << "Gap starts at index = " << startIdx << "\n";
	//cout << "Gap width = " << width << "\n";

	return { minVal, startIdx, width };
}

int findBestFittingRect(const tuple<int, int, int>& niche, const vector<vector<int>>& items2, set<int>& packed_items, vector<tuple<int, int>>& packing_locations) {
	//while finding the best rectangle three posibilities:  
	//(i)there is a shape with a width matching the gap - if multiple, choose the one width largest height one
	//(ii) width smaller than the gap
	//(iii) no item
	//As soon as a match is found, terminate

	int minHeight = get<0>(niche);
	int startIdx = get<1>(niche);
	int gapWidth = get<2>(niche);

	int packed_item_index = -1;
	for (int i = 0; i < items2.size(); i++) {
		if (packed_items.count(i)) continue; //i is already packed

		if (items2[i][0] <= gapWidth) {
			//cout << "\n Found an item!" << endl;

			//use left-most placement policy
			packing_locations[i] = { startIdx , minHeight }; //update the packing coordinate (bottom left)
			//cout << "Item-" << i << " is placed at (" << startIdx << "," << minHeight << ")" << endl;
			packed_items.insert(i); //update packed items set
			packed_item_index = i;
			return packed_item_index;
		}

	}
	return packed_item_index;
}

void BEST_FIT(Instance& inst, Solution& sol) {
	//cout << "Strip W = " << inst.W << endl;
	//Items are already sorted from largest width to smallest (ties: largest height to smallest)
	// from CSP to BPP
	vector<vector<int>> items2;
	int n2 = 0;
	for (int j = 0; j < inst.n; j++) {
		for (int k = 0; k < inst.items[j][2]; k++) {
			items2.push_back({ inst.items[j][0],inst.items[j][1],j });
			//cout << "item=" << n2 << ": " << inst.items[j][0] << ", " << inst.items[j][1] << ", " << j << endl;
			n2++;
		}
	}
	inst.items2 = items2;

	//Create a skyline array to keep track of the height in each column
	vector<int> skyline_array(inst.W, 0);
	set<int> packed_items;
	vector<tuple<int, int>> packing_locations(n2, { -1,-1 }); // initialize all packing locations as (-1,-1)

	while (packed_items.size() < n2) {
		//Find the lowest gap (niche), check the skyline_array, find the smallest valued entry index, and for the size of the gap count the consecutive array items of equal value
		auto niche = findLowestGap(skyline_array);
		int minHeight = get<0>(niche);
		int startIdx = get<1>(niche);
		int gapWidth = get<2>(niche);

		//Find the best-fitting rectangle
		int bestRect = findBestFittingRect(niche, items2, packed_items, packing_locations);

		if (bestRect != -1) {
			//Raise the skyline
			for (int col = startIdx; col < startIdx + items2[bestRect][0]; col++)  skyline_array[col] += items2[bestRect][1];
		}
		else {
			//no item could be packed, wasted space! Raise the skyline to the lowest neighbour
			int leftneighbour = 99999, rightneighbour = 99999;

			if (startIdx != 0) leftneighbour = skyline_array[startIdx - 1];
			if (startIdx + gapWidth != inst.W) rightneighbour = skyline_array[startIdx + gapWidth];

			int chosen_height = 0;
			if (leftneighbour <= rightneighbour) chosen_height = leftneighbour;// raise the gap to the size of left neighbour
			else chosen_height = rightneighbour; // raise the gap to the size of right neighbour
			// raise the gap to the size of chosen height
			//cout << "no item can be packed in the lowest gap! lowest gap raised to " << chosen_height << endl;
			for (int col = startIdx; col < startIdx + gapWidth; col++)  skyline_array[col] = chosen_height;

		}
	}

	sol.UB = *max_element(skyline_array.begin(), skyline_array.end());
	sol.packing_locations = packing_locations;
	//cout << "Heuristic solution = " << sol.UB + sol.prepacked_height << endl;
	return;
}
