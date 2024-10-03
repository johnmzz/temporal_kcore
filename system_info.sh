#!/bin/bash

echo "System Information" > system_info.txt
uname -a >> system_info.txt

echo -e "\nCPU Information" >> system_info.txt
lscpu >> system_info.txt

echo -e "\nMemory Information" >> system_info.txt
free -h >> system_info.txt

echo -e "\nDisk Information" >> system_info.txt
df -h >> system_info.txt
lsblk >> system_info.txt

echo -e "\nNetwork Information" >> system_info.txt
ip a >> system_info.txt

echo -e "\nOperating System Information" >> system_info.txt
cat /etc/os-release >> system_info.txt

echo -e "\nKernel Messages" >> system_info.txt
dmesg | tail -n 50 >> system_info.txt  # Only get the last 50 lines for brevity

echo -e "\nPCI Devices" >> system_info.txt
lspci >> system_info.txt

echo -e "\nUSB Devices" >> system_info.txt
lsusb >> system_info.txt

echo -e "\nRunning Processes" >> system_info.txt
ps aux >> system_info.txt

# View the collected system information
cat system_info.txt

# to run this bash script, in kernel run:
# chmod +x get_system_info.sh
# ./get_system_info.sh