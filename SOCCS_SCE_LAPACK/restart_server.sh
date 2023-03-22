if [ -z "$1" ];then
  echo "Usage: $0 <index_0> <index_1> ..."
  exit 255
fi

killall server_soccs
for index in "$@"
do
   echo "Start server at core $index -"
   ./server_soccs $index &
done

