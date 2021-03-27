# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 1 - dpchange**

This code returns the min number of coins necessary to reach an amount of money given the types of coins available  (non-exhaustive result)

```python
def dpchange(money, coins):
    num_coins = len(coins)
    max_coin = max(coins)
    min_num_coins = {0:0}
    coins_list = {0:[]*num_coins}
    for m in range(1, money+1):
        min_num_coins[m] = money/min(coins) + 1
        if m <= max_coin:
            for coin in coins:
                if m >= coin:
                    if min_num_coins[m-coin] + 1 < min_num_coins[m]:
                        min_num_coins[m] = min_num_coins[m-coin] + 1
                        k = coins_list[m-coin].copy()
                        c = {m:k}
                        coins_list.update(c)
                        coins_list[m].append(coin)

         else:
            min_num_coins.pop(m-max_coin-1)
            coins_list.pop(m-max_coin-1)
            for coin in coins:
                if m >= coin:
                    if min_num_coins[m-coin] + 1 < min_num_coins[m]:
                        min_num_coins[m] = min_num_coins[m-coin] + 1
                        coins_list[m] = coins_list[m-coin]
                        coins_list[m].append(coin)
            
    return min_num_coins[money], coins_list[money]
 ```
