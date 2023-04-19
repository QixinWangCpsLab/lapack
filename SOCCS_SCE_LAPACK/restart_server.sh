if [ -z "$1" ];then
  echo "Usage: $0 <index_0> <index_1> ..."
  exit 255
fi

killall server_soc
make -s
for index in "$@"
do
   echo "Start server at core $index -"
   sudo chrt -f 99 ./server_soc $index &
done

