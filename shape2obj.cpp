#include <iostream>
#include <string>
using namespace std;

int main()
{
	string line;
	getline(cin, line);
	int n = stoi(line.substr(0, line.rfind(" ")));
	for (int i = 0; i < n; i++)
	{
		getline(cin, line);
		cout << "v " << line << endl;
	}
	while (getline(cin, line))
	{
		cout << "f " << line << endl;
	}

	return 0;
}