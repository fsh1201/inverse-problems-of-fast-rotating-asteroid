/* This program calls the program period_scan to look for the best spin state, using multi-threads.

   syntax:
   spin_axis lc input_period_scan out_periods
*/

/*
Copyright (C) 2023  Shuai Feng

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <iostream>
#include <string>
#include <sstream>
#include <thread>

using namespace std;

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		printf("spin_axis lc input_period_scan out_periods\n");
		exit(0);
	}

	int nf;
	thread th[42];
	for (nf = 0; nf < 21; nf++)
	{
		stringstream cmd("");
		cmd << "cat " << string(argv[1]) << " | ./period_scan " << string(argv[2]) << " " << string(argv[3]) << "_" << nf << ".txt -t " << nf;
		cout << cmd.str() << endl;
		th[nf] = thread([&cmd]() {
			system(cmd.str().c_str());
			});
	}
	for (nf = 0; nf < 21; nf++)
	{
		th[nf].join();
	}

	for (nf = 21; nf < 42; nf++)
	{
		stringstream cmd("");
		cmd << "cat " << string(argv[1]) << " | ./period_scan " << string(argv[2]) << " " << string(argv[3]) << "_" << nf << ".txt -t " << nf;
		cout << cmd.str() << endl;
		th[nf] = thread([&cmd]() {
			system(cmd.str().c_str());
			});
	}
	for (nf = 21; nf < 42; nf++)
	{
		th[nf].join();
	}

	stringstream cmd("");
	cmd << "sort -n -k5 " << string(argv[3]) << " > " << string(argv[3]) << "_out_results.txt";
	cout << cmd.str() << endl;
	system(cmd.str().c_str());

	return 0;
}
