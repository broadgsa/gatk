echo "Running"
bjobs -u gsaadm -r | grep gsaadm | awk '{print $4}' | sort | uniq -c

echo "Pending"
bjobs -u gsaadm -p | grep gsaadm | awk '{print $4}' | sort | uniq -c

echo "Suspended"
bjobs -u gsaadm -s | grep gsaadm | awk '{print $4}' | sort | uniq -c
