#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h> 
#include <sys/wait.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		printf("spin_axis lc input_period_scan out_periods\n");
		exit(0);
	}
	int nf;
	char cmd[1024];
	pid_t fpid[42], ppid;
	ppid = getpid();
	for (nf = 0; nf < 21; nf++)
	{
		if (getpid() == ppid)
		{
			fpid[nf] = fork();
			if (fpid[nf] < 0)
			{
				printf("error in fork!");
			}
			else if (fpid[nf] == 0)
			{
				// printf("Child process id is %d\n",getpid());   
				sprintf(cmd, "cat %s | /mnt/d/wuli/fast_P/version_0.2.1_f/convexinv/period_scan %s %s_%d.txt -t %d", argv[1], argv[2], argv[3], nf, nf);
				// printf("%s\n", cmd);
				system(cmd);
				exit(0);
			}
			else
			{
				// printf("Parent process id is %d\n",getpid());   
			}
		}
	}
	for (nf = 0; nf < 21; nf++)
	{
		int status;
		waitpid(fpid[nf], &status, 0);
	}
	ppid = getpid();
	for (nf = 21; nf < 42; nf++)
	{
		if (getpid() == ppid)
		{
			fpid[nf] = fork();
			if (fpid[nf] < 0)
			{
				printf("error in fork!");
			}
			else if (fpid[nf] == 0)
			{
				// printf("Child process id is %d\n",getpid());   
				sprintf(cmd, "cat %s | /mnt/d/wuli/fast_P/version_0.2.1_f/convexinv/period_scan %s %s_%d.txt -t %d", argv[1], argv[2], argv[3], nf, nf);
				// printf("%s\n", cmd);
				system(cmd);
				exit(0);
			}
			else
			{
				// printf("Parent process id is %d\n",getpid());   
			}
		}
	}
	for (nf = 21; nf < 42; nf++)
	{
		int status;
		waitpid(fpid[nf], &status, 0);
	}
	for (nf = 0; nf < 42; nf++)
	{
		memset(cmd, 0, 1024);
		sprintf(cmd, "cat %s_%d.txt >> %s", argv[3], nf, argv[3]);
		system(cmd);
		memset(cmd, 0, 1024);
		sprintf(cmd, "rm %s_%d.txt", argv[3], nf);
		system(cmd);
	}
	memset(cmd, 0, 1024);
	sprintf(cmd, "sort -n -k5 %s > %s_out_results.txt", argv[3], argv[3]);
	system(cmd);

	return 0;
}